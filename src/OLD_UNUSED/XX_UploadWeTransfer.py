# conda create -n wetransfer python=2.7
# conda activate wetransfer
# conda install pip
# pip install wetransferpy
# pip install wetransfer-upload
# 

# source ~/code/tfcf/setup.sh

#from wetransferpy import WeTransfer
from wetransfer import WeTransfer

import os

wt = WeTransfer()
dataPath = os.getenv('HOME')
wt.start(dataPath + "/test.txt")

dataPath = os.getenv('DATA')
dataPath = os.getenv('HOME')

wt = WeTransfer(username="hodpa@haem.cam.ac.uk",
                password="Huntly22",
                sender="nikolaus.fortelny@plus.ac.at",
                receivers=["nikolausfortelny@gmail.com"],
                channel='',
                message='Test1',
                expire_in='3m',
                progress=True,
)
wt.uploadFile(dataPath + "/test.txt")

os.listdir(dataPath)