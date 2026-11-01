@echo off
rem Serves the web app (web/) on http://127.0.0.1:8100 and opens it in the default browser.
cd /d "%~dp0"
start "" http://127.0.0.1:8100/
node web\serve.mjs 8100
