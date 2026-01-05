@echo off
setlocal

rem ============================================
rem 引数:
rem   %1 = INP のフルパス
rem   %2 = NCPUS（並列コア数, 例: 4）
rem   %3 = LOGDIR（steps.dat を置くフォルダ)
rem   %4 = gms_progress.ps1 のフルパス
rem 出力:
rem   .out は INP と同じフォルダ
rem   *.steps.dat は LOGDIR に <job>.steps.dat
rem ============================================

set "INP=%~1"
set "NCPUS=%~2"
set "LOGDIR=%~3"
set "PS1=%~4"

if "%INP%"=="" (
  echo ERROR: INP fullpath is not specified.
  goto :END
)

if "%NCPUS%"=="" set "NCPUS=8"
if "%LOGDIR%"=="" set "LOGDIR=%CD%\logs"

rem INP の親フォルダ/ベース名
set "INPDIR=%~dp1"
set "JOB=%~n1"

rem 出力パス
set "LOG=%INPDIR%%JOB%.out"
set "STEPLOG=%LOGDIR%\%JOB%.steps.dat"

rem 入力・スクリプト存在確認
if not exist "%INP%" (
  echo ERROR: INP not found: "%INP%"
  goto :END
)

if not exist "%PS1%" (
  echo ERROR: PS1 not found: "%PS1%"
  goto :END
)

rem ログ削除・ディレクトリ作成
if exist "%LOG%" del "%LOG%"
if not exist "%LOGDIR%" mkdir "%LOGDIR%"
if exist "%STEPLOG%" del "%STEPLOG%"

rem GAMESS実行
set "GMSROOT=C:\Users\Public\gamess-64"
set "VERSION=2023.R1.intel"

echo.
echo === RUNNING GAMESS ===
echo INP  : "%INP%"
echo LOG  : "%LOG%"
echo NCPUS: %NCPUS%
echo.

pushd "%GMSROOT%"
rem 標準出力を LOG にリダイレクトして同期実行
call rungms "%INP%" %VERSION% %NCPUS% > "%LOG%" 2>&1
set "GMS_ERR=%ERRORLEVEL%"
popd

echo.
echo GAMESS finished with ERRORLEVEL=%GMS_ERR%
echo OUT:   "%LOG%"


rem gms_progress.ps1 を同期実行

echo.
echo === POST-PROCESSING WITH gms_progress.ps1 ===
echo PS1 : "%PS1%"
echo STEP: "%STEPLOG%"
echo.

powershell -ExecutionPolicy Bypass -File "%PS1%" "%LOG%" "%STEPLOG%"
set "PS_ERR=%ERRORLEVEL%"

echo.
echo gms_progress.ps1 finished with ERRORLEVEL=%PS_ERR%
echo STEPS: "%STEPLOG%"

:END
endlocal & exit /b %GMS_ERR%
