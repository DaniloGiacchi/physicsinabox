from visual.graph import *
import numpy
import sys

# Settings
fullscreen = False
redcyan_glasses = False

show_forces = False
show_velocities = False
show_histograms = True


# Init the window
scene.x = 0
scene.y = 0
scene.width = 1300
scene.height = 1050
scene.range = 1.5
scene.up = (0, 0, 1)
scene.forward = (1, 0, 0)
scene.title = 'physicsinabox'
scene.fullscreen = fullscreen
if redcyan_glasses:
    scene.stereo = 'redcyan'

global paused
paused = False
paused_label = label(text='Paused', visible=False)


num_particles = 80
dt = 0.002
max_force = 100

A = 1
B = 300000

box_length = 1.0


positions = (numpy.random.rand(num_particles, 3) - 0.5)
# positions[0] = [0, 0, 0]
#positions = numpy.array([[0, 0, 0], [0.05, 0.05, 0.05]])

velocities = (numpy.random.rand(num_particles, 3) - 0.5) * 5
#velocities = numpy.zeros((num_particles, 3))
# velocities[0] = [0, 0, 0]
forces = numpy.zeros((num_particles, 3))

#radii = numpy.random.rand(num_particles, 1) / 100
radii = numpy.zeros((num_particles, 1))
#radii[:] = 0.01
radii[:num_particles/2] = 0.01
radii[num_particles/2:] = 0.02
# radii = numpy.array([[0.08], [0.003], [0.01], [0.006], [0.012]])
mass_density = 100000.
masses = 4. / 3. * numpy.pi * radii**3 * mass_density
charges = numpy.zeros((num_particles, 1))
#charges[:num_particles/2] = 1.
#charges[num_particles/2:num_particles/2*2] = -1.
#charges[:3] = 0.5
#charges[3:] = -0.5

connect_points_closer_than = 0.#15


surrounding_box = box(pos=(0, 0, 0), length=box_length, width=box_length, height=box_length, opacity=0.2)


force_vectors = []
velocity_vectors = []
balls = []
for i in range(num_particles):
    force_vectors.append(arrow(pos=positions[i], axis=(0, 0, 0), shaftwidth=0.003, color=color.blue, fixedwidth=True))
    velocity_vectors.append(arrow(pos=positions[i], axis=(0, 0, 0), shaftwidth=0.003, color=color.yellow, fixedwidth=True))
    if charges[i] > 0:
        ball_color = color.red
    elif charges[i] < 0:
        ball_color = color.blue
    else:
        ball_color = color.yellow
    balls.append(sphere(pos=positions[i], color=ball_color, radius=radii[i]))

if show_histograms:
    plot_heights = 250

    gdisplay(x=1300, y=0, height=plot_heights, title='Distance from Center', xtitle='d', ytitle='count')
    position_histogram = ghistogram(bins=numpy.linspace(0, 0.75, 20), color=color.red, accumulate=False)

    gdisplay(x=1300, y=plot_heights, height=plot_heights, title='Velocity', xtitle='v', ytitle='count')
    velocity_histogram = ghistogram(bins=numpy.linspace(0, 20, 20), color=color.yellow, accumulate=True)

    gdisplay(x=1300, y=2*plot_heights, height=plot_heights, title='Velocity in x direction', xtitle='vx', ytitle='count')
    velocity_x_histogram = ghistogram(bins=numpy.linspace(0, 10, 20), color=color.yellow, accumulate=True)

    gdisplay(x=1300, y=3*plot_heights, height=plot_heights, title='Distance between Particles', xtitle='s', ytitle='count')
    distance_histogram = ghistogram(bins=numpy.linspace(0, 2, 100), color=color.green, accumulate=False)


if connect_points_closer_than:
    new_connections = []
    old_connections = []


while True:
    rate(100)

    # process key events
    while scene.kb.keys:
        key = scene.kb.getkey()
        if key == ' ':
            # global paused
            paused = not paused
            paused_label.visible = paused
        elif key == 'b':
            surrounding_box.visible = not surrounding_box.visible
        else:
            print key


    if not paused:

        for calulation_step in range(1):
            forces[:] = 0


            # ---------------
            # external forces
            # ---------------

            # central force
            # forces += -0.05 * positions / numpy.expand_dims(numpy.clip(numpy.linalg.norm(positions, axis=1), 0.1, sys.maxint)**3, 1)

            # earth gravity
            # TODO: Optimize!
            # f = - 100 * masses * (positions + 0.5)
            # f[:, 0] = 0
            # f[:, 1] = 0
            # forces += f


            # --------
            # friction
            # --------
            #forces += - 3. * velocities

            if show_histograms:
                summed_distance = numpy.zeros(0)

            if connect_points_closer_than:
                old_connections = new_connections
                new_connections = []

            for i in range(1, num_particles):
                difference = positions - numpy.roll(positions, -i, axis=0)
                distance = numpy.expand_dims(numpy.clip(numpy.linalg.norm(difference, axis=1), 0.1, sys.maxint), 1)

                if show_histograms:
                    summed_distance = numpy.concatenate((summed_distance, numpy.linalg.norm(difference, axis=1)))

                if connect_points_closer_than:
                    for j in range(num_particles):
                        if distance[j, 0] <= connect_points_closer_than:
                            new_connections.append(curve(pos=(positions[j], positions[(j+i) % num_particles])))

                # ------------------
                # interaction forces
                # ------------------

                # gravity
                #forces += -0.5 * masses * numpy.roll(masses, -i, axis=0) * difference / distance**3

                # coloumb
                #forces += charges * numpy.roll(charges, -i, axis=0) * difference / distance**3

                # lennard jones
                forces += 1.e-10 * difference * (A / distance**14 - B / distance**8)

            #forces.clip(-1, 1)

            if connect_points_closer_than:
                for connection in old_connections:
                    connection.visible = False
                    del connection

            velocities = velocities + dt * forces / masses
            positions = (positions + dt * velocities)

            # -----------------
            # bounce from walls
            # -----------------
            velocities -= 2 * velocities * ((positions > 0.5) + (positions < -0.5))
            #positions = numpy.clip(positions, -0.5, 0.5)

            if show_histograms:
                distance_histogram.plot(summed_distance)

        # Only update the graphics every few steps.
        for i in range(num_particles):
            balls[i].pos = positions[i]

        if show_forces:
            for i in range(num_particles):
                force_vectors[i].pos = positions[i]
                force_vectors[i].axis = forces[i] / 100
        if show_velocities:
            for i in range(num_particles):
                velocity_vectors[i].pos = positions[i]
                velocity_vectors[i].axis = velocities[i] / 30
        if show_histograms:
            position_histogram.plot(numpy.linalg.norm(positions, axis=1))
            velocity_histogram.plot(numpy.linalg.norm(velocities[:num_particles/2], axis=1)**2)
            velocity_x_histogram.plot(numpy.linalg.norm(velocities[num_particles/2:], axis=1)**2)#numpy.abs(velocities[:, 0]))