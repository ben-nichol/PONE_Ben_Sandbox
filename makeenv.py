import sys
import os
import subprocess

envfile = open("env.sh", "w")

basedir = os.path.dirname(os.path.realpath(__file__))

envfile.write(
    "export PYTHONPATH=${PYTHONPATH}:"
    + basedir
    + ":"
    + basedir
    + "/python_modules"
    + "\n"
)
envfile.write("export PONESRCDIR=" + basedir + "\n")

envfile.close()
