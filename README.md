Installation Notes:

This code has been designed using Eigen, and is tested using cmake for installation. 
While cmake 'should' be a bit clever in finding Eigen, in practice I needed to install
Eigen in the following way:

clone Eigen3 into a directory of your choice:
> git clone https://github.com/eigenteam/eigen-git-mirror.git
> cd eigen-git-mirror
> mkdir build
> cd build
> cmake ..
> make install

After you have Eigen3 installation *should* be as simple as:
> cd EXAMPLE_DIR/src
> ccmake .
> make 

You can contact me at daniel.marchand@gmail.com if there are any issues with installation
