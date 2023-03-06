from os import listdir
from os.path import isfile,isdir,join
from natsort import natsorted
from time import sleep
from scanf import scanf




class ff():
    def __init__(self,p):
        self.sfmt = "NIRPS_R%d_R%d.fits"
        self.lf_fmt = "%s/NIRPS_R%2.2d_R02.fits"
        self.ramp=1
        self.ldir = ''
        self.p = p
        self.wait_for_data()#and update the first ramp
    def find_last_folder(self):
        ls = [f for f in listdir(self.p) if all([f.isnumeric(),isdir(join(self.p,f))])]    
        ls = natsorted(ls)
    
        if len(ls)<1:
            self.ldir = ''
            return            
        self.ldir = join(self.p,ls[-1])
        if len(listdir(self.ldir))==0 and len(ls)>1:
            self.ldir = join(self.p,ls[-2])
           
        print("last folder is %s"%self.ldir)   
        return         
    def update_first_ramp(self):
        '''
        Find the newest ramp number in the last directory
        '''
        ls = [f for f in listdir(self.ldir) if '.fits' in f]
        if len(ls)<1:
            self.ramp = 1
            return            
        ls = natsorted(ls)      
        r = scanf(self.sfmt,ls[-1])[0]
        self.ramp = r
        print("1st ramp is %d"%r)
        return 
    def wait_for_data(self):
        '''1st f() to run. Once it return
        something, we should not call this 
        f() again
        '''
        self.find_last_folder()
        #wait to have a folder to check.
        while(self.ldir==''):
            print("Waithing for data. sleeping 6s")
            sleep(6)
            self.find_last_folder(p)
        
        self.update_first_ramp()
    
    def find_next(self):
        '''
        p: data path, ld: last dir returned by wait_for_data
        '''
        
        #while(1):
        timeout=0
        while(1):
            while(isfile(self.lf_fmt%(self.ldir,self.ramp))==False):
                sleep(6)
                timeout+=1
                if timeout==3:
                    timeout=0
                    break;
            if isfile(self.lf_fmt%(self.ldir,self.ramp)):
                self.ramp+=1
                return join(self.ldir,"NIRPS_R%2.2d_R01.fits"%self.ramp),join(self.ldir,"NIRPS_R%2.2d_R02.fits"%self.ramp)
            self.ramp=1
            self.find_last_folder()
    def find_last(self):
        while(1):
            self.find_last_folder()
            ls = [f for f in listdir(self.ldir) if '.fits' in f]
            ls = natsorted(ls)
            ra,re = scanf(self.sfmt,ls[-1])
            if re==2:
                return join(self.ldir,"NIRPS_R%2.2d_R01.fits"%ra),join(self.ldir,"NIRPS_R%2.2d_R02.fits"%ra)
        
            timeout=0
            while(isfile(join(self.ldir,"NIRPS_R%2.2d_R02.fits"%ra))==0):
                print("waithing for the next read")
                sleep(2)
                timeout+=1;
                if (timeout>3):
                    timeout=0
                    break
            if isfile(join(self.ldir,"NIRPS_R%2.2d_R02.fits"%ra)):
                return join(self.ldir,"NIRPS_R%2.2d_R01.fits"%ra),join(self.ldir,"NIRPS_R%2.2d_R02.fits"%ra)
            
if '__main__' in __name__:      
    

    p = "/opt/HxRG-SERVER/data/"
    ff = ff(p)
    for i in range(10):
        print(ff.find_next())