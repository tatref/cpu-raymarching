#![allow(unused_variables)]

extern crate image;
extern crate cgmath;
extern crate nalgebra as na;


use std::path::Path;
use std::fs::File;

use cgmath::prelude::*;
use cgmath::Vector3;


const MAX_STEPS: usize = 64;


#[derive(Copy, Clone, Debug)]
struct Ray {
    org: Vector3<f32>,
    dir: Vector3<f32>,
}


#[derive(Copy, Clone, Debug)]
struct Camera {
    org: Vector3<f32>,
    dir: Vector3<f32>,
    up: Vector3<f32>,
    fov: f32,
    res: (usize, usize),
}

impl Camera {
    fn get_ray(&self, x: usize, y: usize) -> Ray {
        let x = x as f32;
        let y = y as f32;

        let left = - self.up.cross(self.dir);
        let mut dir = self.dir + self.up * (y * 2. / self.res.1 as f32 - 1.)
            + left * (x * 2. / self.res.0 as f32 - 1.);

        dir = dir / dir.distance(Vector3::new(0., 0., 0.));

        Ray {
            org: self.org,
            dir,
        }
    }
}

fn dist(v: &Vector3<f32>) -> f32 {
    (v.distance(Vector3::new(0., 0., 0.)) - 5.).min(v.y)
}

fn march(r: &Ray) -> f32 {
    let mut current_pos = r.org;
    let mut current_d = 0f32;

    for i in 0..MAX_STEPS {
        let d = dist(&current_pos);
        current_d += d;
        current_pos += d * r.dir;

        if d < 0.1 {
            return current_d;
        }
    }


    std::f32::INFINITY
}



fn main() {
    println!("Hello, world!");

    let c = Camera {
        org: Vector3::new(0f32, 10f32, 0f32),
        dir: Vector3::unit_x(),
        up: Vector3::unit_y(),
        fov: 0.5,
        res: (500, 500),
    };

    let mut imgbuf = image::ImageBuffer::new(c.res.0 as u32, c.res.1 as u32);

    for y in 0..c.res.1 {
        for x in 0..c.res.0 {
            let r = c.get_ray(x, y);

            let dist = march(&r);

            if dist == std::f32::INFINITY {
                *imgbuf.get_pixel_mut(x as u32, y as u32) = image::Luma([0 as u8])
            }
            else {
                if dist < 0. {
                    *imgbuf.get_pixel_mut(x as u32, y as u32) = image::Luma([0 as u8])
                }
                else {
                    *imgbuf.get_pixel_mut(x as u32, y as u32) = image::Luma([dist as u8])
                }
            }
            
        }
    }

    image::ImageLuma8(imgbuf).save("out.png").unwrap();
}
