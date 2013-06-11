# This script removes all of the python bytecode (.pyc) files.  To be run before a git push.

find ./ -iname "*.pyc" -exec rm {} \;
