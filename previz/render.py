import subprocess
import os
import threading

def render_file(start, stop):
    subprocess.run(["./previz", str(start), str(stop)])
    return True

if __name__=="__main__":
    threads=[]
    subprocess.run(["make","-B"])
    for i in range(0,4):
        t = threading.Thread(target=render_file, args=(i*600, (i+1)*600)) # usually i*600, (i+1)*600
        t.start()
        threads.append(t)
    
    for thread in threads:
        thread.join()
        