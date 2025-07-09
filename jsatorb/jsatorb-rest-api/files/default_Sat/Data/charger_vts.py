# -*- coding: utf-8 -*-
# Anthony Le Batteux EAE 130723
import os
import shutil
import sys

currentWorkingDirectory = os.getcwd()
vts_file = os.path.join(currentWorkingDirectory, 'default_Sat.vts')
vts_dir  = "/app/jsatorb-rest-api/files"
vts_file_dir  = "/app/jsatorb-rest-api/files/default_Sat"
if os.path.exists(vts_file):
    shutil.copy2(vts_file, vts_file_dir)
    shutil.copy2(vts_file, vts_dir)
