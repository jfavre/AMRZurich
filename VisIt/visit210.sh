if /usr/bin/test -z "$userenv_visit_210"; then

if /usr/bin/test -z "$PATH"
then
PATH=/local/apps/VisIt/2.10/current/linux-x86_64/bin:/local/apps/VisIt/2.10/bin; export PATH
else
PATH=$PATH:/local/apps/VisIt/2.10/current/linux-x86_64/bin:/local/apps/VisIt/2.10/bin; export PATH
fi

if /usr/bin/test -z "$LD_LIBRARY_PATH"
then
LD_LIBRARY_PATH=/usr/lib/nvidia-352:/local/apps/VisIt/2.10/current/linux-x86_64/lib; export LD_LIBRARY_PATH
else
LD_LIBRARY_PATH=/usr/lib/nvidia-352:$LD_LIBRARY_PATH:/local/apps/VisIt/2.10/current/linux-x86_64/lib; export LD_LIBRARY_PATH
fi

#export VISIT_MESA_LIB=/local/apps/VisIt/2.10/current/linux-x86_64/lib/libOSMesa.so
export PYTHONPATH=/local/apps/VisIt/2.10/current/linux-x86_64/lib/site-packages
userenv_visit_210=y; export userenv_visit_210
export VISITHOME=/local/apps/VisIt/2.10
export VISITARCHHOME=/local/apps/VisIt/2.10/current/linux-x86_64
export VISITPLUGININSTPRI=~/.visit/current/linux-x86_64/plugins
fi
