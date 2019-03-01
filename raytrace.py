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

You will provide your implementation to each of the functions below.
We have left pseudocode in comments for all problems but extra credits.
For extra credit, you will need to modify this and `scene.py` files.

Note: The code structure here is only a suggestion.  Feel free to
modify it as you see fit.  Your final image must match _perfectly_
the reference images.  If your final image is not **exactly** the
same, you will lose points.

Hint: Use the image to debug your code.  Store x,y,z values of ray
directions to make that they are pointing in correct direction. Store
x,y,z values of intersection points to make sure they seem reasonable.
Store x,y,z values of normals at intersections to make sure they are
reasonable.  Etc.

Hint: Implement in stages.  Trying to write the entire raytracing
system at once will often introduce multiple errors that will be
very tricky to debug.  Use the critical thinking skill of Developing
Sub-Goals to attack this project!

Remove the `@put_your_code_here` function decorators to remove the
`>>> WARNING <<<` reports when running your code.
'''


class Intersection:
    '''
    Stores information about the point of intersection:
        ray_t: t value along ray to intersection (evaluating ray at t will give pos)
        frame: shading frame of intersection (o: point of intersection, z: normal at intersection)
        mat:   the material of surface that was intersected
    '''
    
    def __init__(self, ray_t:float, frame:Frame, mat:Material):
        self.ray_t = ray_t
        self.frame = frame
        self.mat   = mat


def intersect(scene:Scene, ray:Ray):
    ''' returns shading frame at intersection of ray with scene; otherwise returns None '''
    
    '''
    foreach surface
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
    return closest intersection
    '''

    intersection = None

    for surface in scene.surfaces:
        o = surface.frame.o
        r = surface.radius
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

def irradiance(scene:Scene, ray:Ray, iterations=0):
    ''' computes irradiance (color) from scene along ray (reversed) '''
    
    '''
    get scene intersection
    if not hit, return background
    accumulate color starting with ambient
    foreach light
        compute light response
        compute light direction
        compute material response (BRDF*cos)
        check for shadows and accumulate if needed
    if material has reflections
        create reflection ray
        accumulate reflected light (recursive call) scaled by material reflection
    return accumulated color
    '''

    intersection = intersect(scene, ray)
    if not intersection:
        return scene.background
    
    final_color = Vector((0, 0, 0))
    final_color += scene.ambient * intersection.mat.kd
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
        
    v = -ray.d
    n = intersection.frame.z
    rd = -v + 2 * (v.dot(n)) * n
    r = Ray(intersection.frame.o, rd)
    if iterations < 2:
        final_color += intersection.mat.kr * irradiance(scene, r, 5)
    
    return final_color
    

@timed_call('raytrace') # <= reports how long this function took
def raytrace(scene:Scene):
    ''' computes image of scene using raytracing '''
    
    image = Image(scene.resolution_width, scene.resolution_height)
    
    '''
    if no anti-aliasing
        foreach image row (scene.resolution_height)
            foreach pixel in row (scene.resolution_width)
                compute ray-camera parameters (u,v) for pixel
                compute camera ray
                set pixel to color raytraced with ray
    else
        foreach image row
            foreach pixel in row
                init accumulated color
                foreach sample in y
                    foreach sample in x
                        compute ray-camera parameters
                        computer camera ray
                        accumulate color raytraced with ray
                set pixel to accum color scaled by number of samples
    return rendered image
    '''

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
                u = (col + 0.5) / (scene.resolution_width)
                v = 1 - ((row + 0.5) / (scene.resolution_height))
                q = o + (u - 0.5) * w * x + (v - 0.5) * h * y - (d * z)
                r = Ray.from_segment_no_max(o, q)
                image[col, row] = irradiance(scene, r)

    elif scene.pixel_samples > 1:
        for col in range(scene.resolution_width):
            for row in range(scene.resolution_height):
                color = Vector((0, 0, 0))
                for col2 in range(scene.pixel_samples):
                    for row2 in range(scene.pixel_samples):
                        u = (col + (col2 + 0.5) / scene.pixel_samples) / scene.resolution_width
                        v = 1 - (row + (row2 + 0.5) / scene.pixel_samples) / scene.resolution_height
                        q = o + (u - 0.5) * w * x + (v - 0.5) * h * y - (d * z)
                        r = Ray.from_segment_no_max(o, q)
                        color += irradiance(scene, r)
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
    scene = scene_from_file(scene_filename)
    image = raytrace(scene)
    image.save(image_filename)


print('Done')
print()
