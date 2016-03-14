# how to built the VisIt plugin
I have run this procedure on my Linux desktop
run whatever initialization sh script you have to correctly initialize VisIt.
Mine is in visit210.sh

xml2cmake AMAZE.xml
xml2plugin AMAZE.xml
ccmake .

pushd /local/apps/VisIt/2.10/current/linux-x86_64/lib
ln -s /usr/lib/x86_64-linux-gnu/libdl.so .
ln -s /usr/lib/x86_64-linux-gnu/libhwloc.so .
popd
make

In order to find hdf5.h, your VisIt installation needs to INSTALL the 3-rd party:
VISIT_INSTALL_THIRD_PARTY=ON
