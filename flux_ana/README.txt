This directory makes use of pyroot functionality.

If you know how to enable pyroot functionality, you probably don't need to read this.
In my case, I have multiple version of python installed and the cheater's way to fix this seems to simply use "python2.7" instead of simply "python".
This seems to be dependent on which version of python pyroot was built on. We'll have to see about my other packages...

To ensure that python knows where to look I set some paths:


export PYTHONPATH=/usr/local/Cellar/root/6.12.04_1/lib/root/
export LD_LIBRARY_PATH=/usr/local/Cellar/root/6.12.04_1/lib/root/
export ROOTSYS=/usr/local/Cellar/root/6.12.04_1/lib/root/
