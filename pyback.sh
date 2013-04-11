#!/bin/bash
LIST='ls *.py'
for i in $LIST;
do 
cp $i ~/Desktop/subversion/newby/python_scripts/
echo "Copying ${i}"
done

# run svn
echo "Moving to svn directory" 
cd ~/Desktop/subversion/newby/python_scripts/
echo " ----- SVN Status ----- "
svn status
echo " ----- SVN Update ----- "
svn update

# Ask to add new python scripts
echo "Add new Python files? (y|n) :"
read add
if [ "$add" == "y" ]; then
echo " ----- SVN Add ----- "
svn add ./*py
else
echo "Skipping add"
fi

# Ask to commit changes
echo "Commit changes? (y|n) :"
read com
if [ "$com" == "y" ]; then
echo " ----- SVN Commit ----- "
svn commit
else
echo "Skipping commit"
fi

# return to Python_Enc folder
cd ~/Dropbox/Research/Python_Enc/
echo "Returning to home directory"

# move images to the image folder
echo "Archive Images? (y|n) :"
read arc
if [ "$arc" == "y" ]; then
echo " ----- Archiving Images ----- "
IMAGES='ls *.png'
for i in $IMAGES;
do 
mv $i ~/Dropbox/Research/Python_Enc/image_holder/
echo "Copying ${i} to /image_holder"
done
IMAGES='ls *.ps'
for i in $IMAGES;
do 
mv $i ~/Dropbox/Research/Python_Enc/image_holder/
echo "Copying ${i} to /image_holder"
done
else
echo "Skipping archive"
fi
echo " --- Done --- "

