#!/usr/bin/env bash

# 1. llança la teva app Shiny en segon pla
Rscript -e "shiny::runApp('/Users/jfibla/Library/CloudStorage/Dropbox/GWAS_inspector/app.R', host='0.0.0.0', port=8080)" &
APP_PID=$!

# 2. espera uns segons perquè R aixequi el servidor
sleep 4

# 3. obre túnel ngrok cap al port 8080
ngrok http 8080

# 4. quan tanquis ngrok, mata l'app (opcional)
kill $APP_PID
