@echo off
setlocal
cd /d "%~dp0"

echo [1/4] Updating pip...
python -m pip install --upgrade pip setuptools wheel
if errorlevel 1 goto :failed

echo [2/4] Installing Windows-compatible molecular dependencies...
python -m pip install -r requirements-windows.txt
if errorlevel 1 goto :failed

echo [3/4] Installing DT4DDS core without optional Unix-only analysis dependencies...
python -m pip install dt4dds==1.1.0 --no-deps
if errorlevel 1 goto :failed

echo [4/4] Verifying DT4DDS molecular backend...
python -c "import sys, dt4dds; from rapidfuzz.distance import Levenshtein; print('Python:', sys.executable); print('DT4DDS', dt4dds.__VERSION__, '- molecular backend ready'); print('DT4DDS path:', dt4dds.__file__)"
if errorlevel 1 goto :failed

echo.
echo Installation completed. Start the app with:
echo python -m streamlit run app.py
exit /b 0

:failed
echo.
echo Installation failed. Copy the complete error output for diagnosis.
exit /b 1
