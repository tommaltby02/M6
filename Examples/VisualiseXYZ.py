import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider, Button
from matplotlib.animation import FuncAnimation
import glob

def plot_xyz_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        data = [[float(num) for num in line.split()] for line in lines]

    data = list(zip(*data))
    return data

def update_plot(val, data, sc):
    sc._offsets3d = (data[int(val)][0], data[int(val)][1], data[int(val)][2])

def update_slider(val, slider):
    current_val = slider.val
    max_val = slider.valmax
    if current_val == max_val:
        slider.set_val(slider.valmin)
    else:
        slider.set_val(current_val + 1)



def visualize_xyz_files_with_pause(folder_path):
    file_list = glob.glob(folder_path + "/*.xyz")
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.axis('off')
    ordered_file_list = []
    for i in range(len(file_list)):
        ordered_file_list.append(folder_path + '/config' + str(i) + '.xyz')

    data = [plot_xyz_file(file_path) for file_path in ordered_file_list]

    sc = ax.scatter(data[0][0], data[0][1], data[0][2], s=100)


    ax_slider = plt.axes([0.1, 0.01, 0.8, 0.03])
    slider = Slider(ax_slider, 'Frame', 0, len(data) - 1, valinit=0, valstep=1)

    def update_frame(val):
        update_plot(val, data, sc)
        ax.set_xlim(min(data[val][0]), max(data[val][0]))  # Update the limits for x-axis
        ax.set_ylim(min(data[val][1]), max(data[val][1]))  # Update the limits for y-axis
        ax.set_zlim(min(data[val][2]), max(data[val][2]))  # Update the limits for z-axis

    slider.on_changed(update_frame)

    def animate(frame):
        update_slider(frame, slider)
        update_frame(frame)

    ani = FuncAnimation(fig, animate, frames=len(data)-1, interval=1)


    plt.show()

folder_path = "XYZFiles"
visualize_xyz_files_with_pause(folder_path)

