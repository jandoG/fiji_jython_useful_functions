#@ImagePlus imp0

from ij import IJ, WindowManager
from ij import ImagePlus
from ij import ImageStack
from ij.plugin import Duplicator, ChannelSplitter
import os 
from os import path 


imp = imp0.duplicate()

imps = ChannelSplitter().split(imp)

impch1 = imps[0]
impch2 = imps[1]

impch1_name = impch1.getTitle()
impch2_name = impch2.getTitle()

impch1.setT(1);
impch1.show()
impch2.show()

IJ.run("MultiStackReg", "stack_1="+ impch1_name +" action_1=Align file_1=[] stack_2="+ impch2_name +" action_2=[Align to First Stack] file_2=[] transformation=Translation");

IJ.run("Merge Channels...", "c1="+ impch1_name +" c2="+ impch2_name +" create");
composite = IJ.getImage();
composite.setTitle("Final-drift corrected")
composite.show()
