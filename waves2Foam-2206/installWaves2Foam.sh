# install necessary packages
sudo apt-get install libgsl-dev
sudo apt install subversion

# download OceanWave3D-Fortran90 and move it to ThirdParty
git clone https://github.com/boTerpPaulsen/OceanWave3D-Fortran90.git 
cp -r OceanWave3D-Fortran90 ThirdParty/OceanWave3D-Fortran90BK

# download waves2Foam
svn co http://svn.code.sf.net/p/openfoam-extend/svn/trunk/Breeder_1.6/other/waves2Foam

# go to waves2Foam and copy my Thirdparty.
cd waves2Foam
rm -rf ThirdParty
cp -f ../ThirdParty .
cd ..

# move waves2Foam to $WM_PROJECT_USER_DIR/applications/utilities/
mkdir -p $WM_PROJECT_USER_DIR/applications/utilities/
cp -r waves2Foam $WM_PROJECT_USER_DIR/applications/utilities/

# compile
cd $WM_PROJECT_USER_DIR/applications/utilities/waves2Foam

./Allwmake
