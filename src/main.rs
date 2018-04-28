#![allow(unused_imports)]
#![allow(unused_variables)]

extern crate image;
extern crate num_traits;
extern crate nalgebra as na;


use std::f32;
use std::path::Path;
use std::fs::File;

use num_traits::Pow;
use na::{Vector3, Point3, Unit};
use num_traits::identities::Zero;


const MAX_STEPS: usize = 64;
const COLLISION_DISTANCE: f32 = 0.01;



trait Sdf {
    fn sdf(&self, v: &Vector3<f32>) -> f32;
}


#[derive(Copy, Clone, Debug)]
struct Ray {
    org: Vector3<f32>,
    dir: Unit<Vector3<f32>>,
}


#[derive(Copy, Clone, Debug)]
struct Camera {
    org: Vector3<f32>,
    dir: Unit<Vector3<f32>>,
    up: Unit<Vector3<f32>>,
    fov: f32,
    res: (usize, usize),
}

impl Camera {
    fn get_ray(&self, x: usize, y: usize) -> Ray {
        let x = x as f32;
        let y = y as f32;

        let left = - self.up.cross(&self.dir);
        let dir = *self.dir + *self.up * (y * 2. / self.res.1 as f32 - 1.)
            + left * (x * 2. / self.res.0 as f32 - 1.);


        Ray {
            org: self.org,
            dir: Unit::new_normalize(dir),
        }
    }
}

struct Sphere {
    center: Vector3<f32>,
    radius: f32,
}

impl Sdf for Sphere {
    fn sdf(&self, v: &Vector3<f32>) -> f32 {
        (v - self.center).norm() - self.radius
    }
}

struct Plane {
    normal: Unit<Vector3<f32>>,
    dist: f32,
}

impl Sdf for Plane {
    fn sdf(&self, v: &Vector3<f32>) -> f32 {
        v.dot(&*self.normal) + self.dist
    }
}


struct Scene {
    camera: Camera,
    sdf: Fn(&Vector3<f32>) -> f32,
}

impl Scene {
    fn render(&self) {
    }
}

fn dist(v: &Vector3<f32>) -> f32 {
    let sphere = Sphere { center: Vector3::new(10., 0., 0.), radius: 5f32 };

    let plane = Plane { normal: Unit::new_normalize(Vector3::new(0., 1., 0.)), dist: 0. };

    plane.sdf(v).min(sphere.sdf(v))
}

fn march(r: &Ray) -> f32 {
    let mut current_pos = r.org;
    let mut current_d = 0f32;

    for i in 0..MAX_STEPS {
        let d = dist(&current_pos);

        if d < COLLISION_DISTANCE {
            return current_d;
        }

        current_d += d;
        current_pos += *r.dir * d;
    }

    std::f32::INFINITY
}



fn main() {
    println!("Hello, world!");

    let c = Camera {
        org: Vector3::new(-10., 10., 0.),
        dir: Vector3::x_axis(),
        up: Vector3::y_axis(),
        fov: 0.5,
        res: (500, 500),
    };

    let back_color: image::Rgb<u8> = image::Rgb([200, 200, 200]);

    let mut imgbuf = image::ImageBuffer::new(c.res.0 as u32, c.res.1 as u32);

    let clock = std::time::Instant::now();

    for y in 0..c.res.1 {

        for x in 0..c.res.0 {
            let r = c.get_ray(x, y);
            let y = c.res.1 - y -1;  // flip vertficaly

            let dist = march(&r);
            let pos = r.org + *r.dir * dist;

            if dist == std::f32::INFINITY {
                *imgbuf.get_pixel_mut(x as u32, y as u32) = back_color;
            }
            else {
                if dist < 0. {
                    *imgbuf.get_pixel_mut(x as u32, y as u32) = image::Rgb([0, 0, 0]);
                }
                else {
                    let color = [
                        (pos.x.fract() * 256f32) as u8,
                        (pos.y.fract() * 256f32) as u8,
                        (pos.z.fract() * 256f32) as u8,
                        
                    ];

                    *imgbuf.get_pixel_mut(x as u32, y as u32) = image::Rgb(color);
                }
            }

            if x < 10 && y < 10 {
                *imgbuf.get_pixel_mut(x as u32, y as u32) = image::Rgb([250, 250, 250]);
            }
            
        }
    }
    let elapsed = clock.elapsed().subsec_millis();
    println!("render time: {:?}", elapsed);

    image::ImageRgb8(imgbuf).save("out.png").unwrap();
}
