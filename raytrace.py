import os
import sys
import math

from common.utils import put_your_code_here, timed_call
from common.maths import Vector, Point, Normal, Ray, Direction, Frame, sqrt
from common.scene import Scene, Material, scene_from_file
from common.image import Image

'''
The following functions provide algorithms for raytracing a scene.

`raytrace` renders an image for given scene by calling `irradiance`.
`irradiance` computes irradiance from scene along ray (reversed).
`intersect` computes intersection details between ray and scene.
'''


'''
Stores information about the point of intersection:
    ray_t: t value along ray to intersection (evaluating ray at t will give pos)
    frame: shading frame of intersection (o: point of intersection, z: normal at intersection)
    mat:   the material of surface that was intersected
'''
class Intersection:
    def __init__(self, ray_t:float, frame:Frame, mat:Material):
        self.ray_t = ray_t
        self.frame = frame
        self.mat   = mat


# returns shading frame at intersection of ray with scene; otherwise returns None
def intersect(scene:Scene, ray:Ray):
    intersection = None
    for surface in scene.surfaces:
        o = surface.frame.o
        r = surface.radius
        '''
        if surface is a quad
            compute ray intersection (and ray_t), continue if not hit
            check if computed ray_t is between min and max, continue if not
            check if this is closest intersection, continue if not
            record hit information
        if surface is a sphere
            compute ray intersection (and ray_t), continue if not hit
            check if computed ray_t is between min and max, continue if not
            check if this is closest intersection, continue if not
            record hit information
        '''
        if surface.is_quad:
            n = surface.frame.z
            t = (o - ray.e).dot(n) / ray.d.dot(n)
            p = ray.eval(t)
            p_ = surface.frame.w2l_point(p)
            if abs(p_.x) > r or abs(p_.y) > r:
                continue
        elif surface.is_circle:
            n = surface.frame.z
            t = (o - ray.e).dot(n) / ray.d.dot(n)
            p = ray.eval(t)
            p_ = surface.frame.w2l_point(p)
            if (p_.x ** 2 + p_.y ** 2) > r ** 2:
                continue
        else:
            b = 2 * ray.d.dot(ray.e - o)
            c = (ray.e - o).length_squared - r * r
            d = b * b - 4 * c
            if d < 0:
                continue
            t = (-b - math.sqrt(d)) / 2

        if not ray.valid_t(t):
            continue

        if intersection is None or t <= intersection.ray_t:
            p = ray.eval(t)
            n = Ray.from_segment_no_max(o, p).d
            if surface.is_quad:
                n = surface.frame.z
            f = Frame(o=p, z=n)
            m = surface.material
            intersection = Intersection(t, f, m)

    return intersection

# computes irradiance (color) from scene along ray (reversed)
def irradiance(scene:Scene, ray:Ray, iterations=0):
    # get scene intersection
    intersection = intersect(scene, ray)
    # if not hit, return background
    if not intersection:
        return scene.background

    # accumulate color starting with ambient
    final_color = Vector((0, 0, 0))
    final_color += scene.ambient * intersection.mat.kd

    '''
    foreach light
        compute light response
        compute light direction
        compute material response (BRDF*cos)
        check for shadows and accumulate if needed
    '''
    for light in scene.lights:
        s = light.frame.o
        p = intersection.frame.o
        n = intersection.frame.z
        ray_to_light = Ray.from_segment(p, s)
        if intersect(scene, ray_to_light):
            continue
        if light.is_point:
            response = light.intensity / (s - p).length_squared
            direction = (s - p) / (s - p).length
            h = (ray_to_light.d + -ray.d)
            h /= h.length
            final_color += response * (intersection.mat.kd + intersection.mat.ks * max(0, n.dot(h)) ** intersection.mat.n) * max(0, (n.dot(direction)))
        else:
            response = light.intensity
            direction = light.frame.z
            h = (ray_to_light.d + -ray.d)
            h /= h.length
            final_color += response * (intersection.mat.kd + intersection.mat.ks * max(0, n.dot(h)) ** intersection.mat.n) * max(0, (n.dot(direction)))

    # create reflection ray
    # accumulate reflected light (recursive call) scaled by material reflection
    v = -ray.d
    n = intersection.frame.z
    rd = -v + 2 * (v.dot(n)) * n
    r = Ray(intersection.frame.o, rd)
    if iterations < 2:
        final_color += intersection.mat.kr * irradiance(scene, r, 5)

    return final_color


# computes image of scene using raytracing
@timed_call('raytrace') # <= reports how long this function took
def raytrace(scene:Scene):
    image = Image(scene.resolution_width, scene.resolution_height)

    o = scene.camera.frame.o
    x = scene.camera.frame.x
    y = scene.camera.frame.y
    z = scene.camera.frame.z
    w = scene.camera.width
    h = scene.camera.height
    d = scene.camera.dist

    if scene.pixel_samples == 1:
        for row in range(scene.resolution_height):
            for col in range(scene.resolution_width):
                # compute ray-camera parameters (u,v) for pixel
                u = (col + 0.5) / (scene.resolution_width)
                v = 1 - ((row + 0.5) / (scene.resolution_height))
                q = o + (u - 0.5) * w * x + (v - 0.5) * h * y - (d * z)
                # compute camera ray
                r = Ray.from_segment_no_max(o, q)
                # set pixel to color raytraced with ray
                image[col, row] = irradiance(scene, r)

    elif scene.pixel_samples > 1:
        for col in range(scene.resolution_width):
            for row in range(scene.resolution_height):
                # init accumulated color
                color = Vector((0, 0, 0))
                for col2 in range(scene.pixel_samples):
                    for row2 in range(scene.pixel_samples):
                        # compute ray-camera parameters
                        u = (col + (col2 + 0.5) / scene.pixel_samples) / scene.resolution_width
                        v = 1 - (row + (row2 + 0.5) / scene.pixel_samples) / scene.resolution_height
                        q = o + (u - 0.5) * w * x + (v - 0.5) * h * y - (d * z)
                        # computer camera ray
                        r = Ray.from_segment_no_max(o, q)
                        # accumulate color raytraced with ray
                        color += irradiance(scene, r)
                # set pixel to accum color scaled by number of samples
                image[col, row] = color / (scene.pixel_samples ** 2)

    return image


if len(sys.argv) < 2:
    print('Usage: %s <path/to/scenefile0.json> [path/to/scenefile1.json] [...]' % sys.argv[0])
    sys.exit(1)


for scene_filename in sys.argv[1:]:
    base,_ = os.path.splitext(scene_filename)
    image_filename = '%s.png' % base

    print('Reading scene: %s' % scene_filename)
    print('Writing image: %s' % image_filename)

    print('Raytracing...')
    scene_filename = scene_filename.strip()
    scene_filename = scene_filename.replace('..', '')
    scene_filename = scene_filename.replace('/', os.path.sep).replace('\\', os.path.sep)
    scene_filename = os.path.join(os.getcwd(), scene_filename)
    scene = scene_from_file(scene_filename)
    image = raytrace(scene)
    image.save(image_filename)


print('Done')
print()
