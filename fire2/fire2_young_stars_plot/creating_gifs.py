from PIL import Image
#import glob  #use it if you want to read all certain file type
    
imgs=[]
for i in range(596,691): #startig and ending value+1 of the index that identifies different file names or imgs = glob.glob("*.png") can do done as well
    filename=str(i)+"test.png"
    imgs.append(filename)
    print("Scanned image from the File: ",filename,end="\r",flush=True)    

frames = []
for i in imgs:
    new_frame = Image.open(i)
    frames.append(new_frame)
# Save into a GIF file that loops forever
gifname="all_clusters_star_formation.gif"
frames[0].save(gifname, format='GIF',
               append_images=frames[1:],
               save_all=True,
               duration=400, loop=0)    #duration =200 looked good
print("\nSaved the GIF into the File: ",gifname)