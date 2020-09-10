you must compile ParaView from source. Once done, 

cd build

cmake ..

make

export PV_PLUGIN_PATH=\`pwd\`/lib/paraview-5.8/plugins/pvAMAZEReader

and you're set to go.
