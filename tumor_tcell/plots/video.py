import os
import shutil

import cv2
import numpy as np
import matplotlib.pyplot as plt

from tumor_tcell.plots.snapshots import (
    make_snapshots_figure,
    format_snapshot_data,
    get_field_range,
    get_agent_colors,
)

PLOT_WIDTH = 7


def make_snapshot_function(
        data,
        bounds,
        agent_colors=None,
        **kwargs):
    agent_colors = agent_colors or {}
    multibody_agents, multibody_fields = format_snapshot_data(data)

    # make the snapshot plot function
    time_vec = list(multibody_agents.keys())

    # get fields and agent colors
    multibody_field_range = get_field_range(multibody_fields, time_vec)
    multibody_agent_colors = get_agent_colors(multibody_agents)
    multibody_agent_colors.update(agent_colors)

    def plot_single_snapshot(t_index):
        time_indices = np.array([t_index])
        snapshot_time = [time_vec[t_index]]
        fig = make_snapshots_figure(
            time_indices=time_indices,
            snapshot_times=snapshot_time,
            agents=multibody_agents,
            agent_colors=multibody_agent_colors,
            fields=multibody_fields,
            field_range=multibody_field_range,
            n_snapshots=1,
            bounds=bounds,
            default_font_size=12,
            plot_width=PLOT_WIDTH,
            # scale_bar_length=0,
            **kwargs)
        return fig

    return plot_single_snapshot, time_vec


def video_from_images(img_paths, out_file):
    # make the video
    img_array = []
    size = None
    for img_file in img_paths:
        img = cv2.imread(img_file)
        height, width, layers = img.shape
        size = (width, height)
        img_array.append(img)

    out = cv2.VideoWriter(out_file, cv2.VideoWriter_fourcc(*'mp4v'), 15, size)
    for i in range(len(img_array)):
        out.write(img_array[i])
    out.release()


def make_video(
        data,
        bounds,
        step=1,
        agent_colors=None,
        out_dir='out',
        filename='snapshot_vid',
        **kwargs
):
    """Make a video with snapshots across time """

    # make images directory, remove if existing
    out_file = os.path.join(out_dir, f'{filename}.mp4')
    images_dir = os.path.join(out_dir, f'_images')
    if os.path.isdir(images_dir):
        shutil.rmtree(images_dir)
    os.makedirs(images_dir)

    # get the single snapshots function
    snapshot_fun, time_vec = make_snapshot_function(
        data,
        bounds,
        agent_colors=agent_colors,
        **kwargs)

    # make the individual snapshot figures
    img_paths = []
    for t_index in range(0, len(time_vec) - 1, step):
        fig_path = os.path.join(images_dir, f"img{t_index}.jpg")
        img_paths.append(fig_path)

        fig = snapshot_fun(t_index)
        fig.savefig(fig_path, bbox_inches='tight')
        plt.close()

    # make the video
    video_from_images(img_paths, out_file)

    # delete image folder
    shutil.rmtree(images_dir)
