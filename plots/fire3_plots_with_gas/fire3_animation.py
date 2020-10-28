import imageio
images=[]

plot_path="./" #creating a path to store the plots only if it does not exist


snapshot_start=596
snapshot_end=690
for img in range(snapshot_start,snapshot_end+1):
    print(img)
    images.append(imageio.imread(plot_path+"fire3_clusters_snap"+str(img) + '.png'))
    
imageio.mimsave(plot_path+"fire3_animation.gif", images, duration = 1/5)
#########################################
#########################################
####warning !!!!! it caused flickering of colors so it is almost useless for creating gifs when we have many color dots     