#!/bin/bash

bash /app/master_setup.sh

exec streamlit run /app/app.py