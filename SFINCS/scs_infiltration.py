"""Classic SCS-CN infiltration, matching the global SCS sensitivity suite.

GCN250 ARCII -> median reprojection -> S=max(1000/CN-10,0) inches.
Missing CN -> CN100 (S=0, no infiltration); no dry-period recovery.
Binary order follows the existing one-based, column-major sfincs.ind.
No solver invocation and no mutation of hard-linked source inputs.
"""
from pathlib import Path
import hashlib
import json
import os
import re
import tempfile
import warnings

def enabled():
    method=os.environ.get('SFINCS_INFILTRATION','scs_cn').lower()
    if method not in {'scs_cn','none'}:raise ValueError('SFINCS_INFILTRATION must be scs_cn or none')
    return method=='scs_cn'

def settings(path):
    values={}
    for line in path.read_text(encoding='utf-8-sig').splitlines():
        if '=' not in line or line.lstrip().startswith(('#','!')):continue
        key,value=line.split('=',1);key=key.strip().lower()
        if key in values:raise ValueError('Duplicate SFINCS key: '+key)
        values[key]=re.split(r'[#!]',value,1)[0].strip().strip('"\'')
    return values

def cn_source(cn_raster=None):
    default=Path(__file__).resolve().parents[1]/'external/soil/GCN250_ARCII.tif'
    path=Path(cn_raster or os.environ.get('SFINCS_GCN250_PATH',default)).expanduser().resolve()
    if not path.is_file():raise FileNotFoundError('Set SFINCS_GCN250_PATH to GCN250_ARCII.tif: '+str(path))
    return path

def atomic(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fd,name=tempfile.mkstemp(dir=path.parent,prefix='scs_',suffix='.tmp')
    try:
        with os.fdopen(fd,'wb') as stream:stream.write(data)
        os.replace(name,path)
    finally:
        if os.path.exists(name):os.unlink(name)

def apply_scs_infiltration(model_dir,cn_raster=None):
    if not enabled():return None
    import numpy as np
    import rasterio
    from rasterio.transform import Affine
    from rasterio.warp import reproject,Resampling
    model=Path(model_dir).resolve();inp=model/'sfincs.inp'
    if (model/'sfincs_map.nc').exists():raise RuntimeError('Refusing to change a simulated model: '+str(model))
    cn=cn_source(cn_raster)
    v=settings(inp)
    for key in ['qinffile','netqinf_file','ksfile','psifile','sigmafile','f0file','fcfile','kdfile','smaxfile','sefffile','inffile','inftype']:
        if v.get(key,'none').lower() not in {'none','nil','','0','-999'}:raise ValueError('Conflicting infiltration setting: '+key)
    if float(v.get('qinf','0'))!=0:raise ValueError('Cannot combine SCS with nonzero qinf')
    keys=['nmax','mmax','dx','dy','x0','y0','rotation','epsg']
    grid={k:v[k] for k in keys};nx,ny=int(v['mmax']),int(v['nmax'])
    index=(model/v['indexfile']).resolve();raw=index.read_bytes()
    idx=np.frombuffer(raw,dtype='<u4');count=int(idx[0]);ids=idx[1:].astype('int64')-1
    if len(ids)!=count or not count or ids.min()<0 or ids.max()>=nx*ny or len(np.unique(ids))!=count:raise ValueError('Invalid sfincs.ind')
    signature=dict(version=1,grid=grid,index_sha256=hashlib.sha256(raw).hexdigest(),cn_path=str(cn),cn_bytes=cn.stat().st_size,cn_mtime_ns=cn.stat().st_mtime_ns)
    key=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()
    cache=Path(os.environ.get('SFINCS_SCS_CACHE',Path(__file__).resolve().parents[1]/'work/cache/scs'))/key
    target=cache/'sfincs.scs';meta=cache/'metadata.json'
    prior=json.loads(meta.read_text()) if meta.exists() and target.exists() else None
    if prior and prior['signature']==signature and hashlib.sha256(target.read_bytes()).hexdigest()==prior['sha256']:
        payload=target.read_bytes();report=prior
    else:
        angle=np.deg2rad(float(v['rotation']));dx,dy=float(v['dx']),float(v['dy'])
        transform=Affine(dx*np.cos(angle),-dy*np.sin(angle),float(v['x0']),dx*np.sin(angle),dy*np.cos(angle),float(v['y0']))
        raster=np.full((ny,nx),np.nan,dtype='float32')
        with rasterio.Env(GDAL_CACHEMAX=128):
            with rasterio.open(cn) as src:
                if src.nodata!=255:raise ValueError('Expected GCN250 nodata=255')
                reproject(rasterio.band(src,1),raster,src_transform=src.transform,src_crs=src.crs,src_nodata=255,dst_transform=transform,dst_crs='EPSG:'+v['epsg'],dst_nodata=np.nan,resampling=Resampling.med,num_threads=1)
        values=raster.T.ravel()[ids];missing=~np.isfinite(values)
        if not np.all((values[~missing]>=1)&(values[~missing]<=100)):raise ValueError('CN must be 1..100')
        storage=np.maximum(1000./np.where(missing,100,values)-10,0).round(3).astype('<f4')
        payload=storage.tobytes()
        report=dict(signature=signature,sha256=hashlib.sha256(payload).hexdigest(),active_cells=count,missing_CN_cells=int(missing.sum()),method='GCN250 ARCII; median; S inches rounded 3 decimals; missing CN100',sfacinf=.2,storecumprcp=1)
        atomic(target,payload);atomic(meta,json.dumps(report,indent=2).encode())
    if len(payload)!=count*4:raise ValueError('Wrong SCS binary length')
    if report['missing_CN_cells']:warnings.warn(f"{report['missing_CN_cells']}/{count} active cells lack CN; S=0 (no infiltration)",RuntimeWarning)
    atomic(model/'sfincs.scs',payload)
    # Replace rather than edit in place: sfincs.inp may be an NTFS/POSIX hardlink.
    remove={'scsfile','sfacinf','storecumprcp','qinf'}
    lines=[line for line in inp.read_text(encoding='utf-8-sig').splitlines() if not ('=' in line and line.split('=',1)[0].strip().lower() in remove)]
    lines+=['scsfile = sfincs.scs','sfacinf = 0.2','storecumprcp = 1']
    atomic(inp,('\n'.join(lines)+'\n').encode())
    atomic(model/'infiltration_metadata.json',json.dumps(report,indent=2).encode())
    return report
