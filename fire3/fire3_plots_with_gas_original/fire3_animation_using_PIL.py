from PIL import Image
import glob  #use it if you want to read all certain file type
plot_path="./plots_pkl/animation_plots/" 

    
imgs=[]
for i in range(596,691): #startig and ending value+1 of the index that identifies different file names or imgs = glob.glob("*.png") can do done as well
    imgs.append("fire3_clusters_snap"+str(i)+'.png')
    print("scanned image",i)    

frames = []
for i in imgs:
    new_frame = Image.open(i)
    frames.append(new_frame)
# Save into a GIF file that loops forever
frames[0].save('fire3_PIL.gif', format='GIF',
               append_images=frames[1:],
               save_all=True,
               duration=200, loop=0)    #duration =200 looked good