@echo off
setlocal
cd /d "%~dp0"
echo TCSPC Workbench 0.1.0
echo Run this from an Anaconda Prompt with your workbench environment activated.
echo.
python -c "import sys; print('Python:', sys.executable); import numpy, scipy.optimize, scipy.signal, matplotlib.pyplot, pandas, streamlit, threadpoolctl; from packaging.version import Version; assert Version(streamlit.__version__) >= Version('1.51'), 'Streamlit 1.51 or newer is required'"
if errorlevel 1 (
    echo.
    echo Setup check failed. Follow QUICK_START_WINDOWS.md or README.md.
    echo Activate your working environment before starting the app.
    pause
    exit /b 1
)
python -m streamlit run app.py
pause
