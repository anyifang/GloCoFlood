"""Run with: python -m unittest discover -s SFINCS -p test_scs_infiltration.py"""
import os
from pathlib import Path
import tempfile
import unittest
import sys
from unittest.mock import patch
import numpy as np
import rasterio
from rasterio.transform import Affine
from scs_infiltration import apply_scs_infiltration,settings

class SCSChecks(unittest.TestCase):
    def test_storage_order_missing_and_hardlink(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);model=root/'model';model.mkdir();cn=root/'cn.tif'
            with rasterio.open(cn,'w',driver='GTiff',height=2,width=2,count=1,dtype='uint8',crs='EPSG:32650',transform=Affine(200,0,500000,0,200,2000000),nodata=255) as d:
                d.write(np.array([[50,100],[75,255]],dtype='uint8'),1)
            np.array([4,1,2,3,4],dtype='<u4').tofile(model/'sfincs.ind')
            original=root/'source.inp'
            original.write_text('nmax=2\nmmax=2\ndx=200\ndy=200\nx0=500000\ny0=2000000\nrotation=0\nepsg=32650\nindexfile=sfincs.ind\nqinf=0\n')
            before=original.read_bytes();os.link(original,model/'sfincs.inp')
            with patch.dict(os.environ,{'SFINCS_INFILTRATION':'scs_cn','SFINCS_SCS_CACHE':str(root/'cache')}):
                report=apply_scs_infiltration(model,cn)
                np.testing.assert_array_equal(np.fromfile(model/'sfincs.scs',dtype='<f4'),np.array([10,3.333,0,0],dtype='<f4'))
                self.assertEqual(report['missing_CN_cells'],1)
                self.assertEqual(original.read_bytes(),before)
                self.assertEqual(settings(model/'sfincs.inp')['storecumprcp'],'1')
                self.assertEqual(apply_scs_infiltration(model,cn)['sha256'],report['sha256'])
                (model/'sfincs_map.nc').touch()
                with self.assertRaises(RuntimeError):apply_scs_infiltration(model,cn)

    def test_disable_is_explicit_noop(self):
        with patch.dict(os.environ,{'SFINCS_INFILTRATION':'none'}):
            self.assertIsNone(apply_scs_infiltration('nonexistent-model'))

    def test_decomposition_keeps_same_soil(self):
        sys.path.insert(0,str(Path(__file__).parent/'global_build/utilities'))
        import build_sfincs_forcing_decomposition_cases as d
        with tempfile.TemporaryDirectory() as temp:
            source=Path(temp)/'source';target=Path(temp)/'no_rain';source.mkdir();target.mkdir()
            (source/'sfincs.scs').write_bytes(b'soil')
            inp=source/'sfincs.inp';inp.write_text('scsfile=sfincs.scs\nsfacinf=0.2\nstorecumprcp=1\nnetamprfile=rain.nc\n')
            result=d.render_inp(d.read_inp(inp),source,target,'no_rain',None)
            self.assertNotIn('netamprfile',result)
            values=dict((line.split('=',1)[0].strip(),line.split('=',1)[1].strip())for line in result.splitlines())
            self.assertEqual((target/values['scsfile']).resolve(),(source/'sfincs.scs').resolve())
            self.assertEqual(values['storecumprcp'],'1')

if __name__=='__main__':unittest.main()
