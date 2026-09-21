@echo off
setlocal
cd /d "%~dp0"
python example\PRD_single_TC\run_prd_single_tc.py %*
exit /b %errorlevel%
