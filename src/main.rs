#![allow(unused_imports)]
#![feature(euclidean_division)]
#![allow(unused_variables)]

extern crate image;
extern crate num_traits;
extern crate nalgebra as na;
extern crate rayon;


use std::f32;
use std::path::Path;
use std::fs::File;

use image::ConvertBuffer;
use na::{Vector3, Vector2, Point3, Unit};
use na::{U1, U2, U3};
use num_traits::identities::Zero;
use num_traits::Pow;
use rayon::prelude::*;



const MAX_STEPS: usize = 256;
const PRECISION: f32 = 0.001;



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
        //(v - self.center).norm() - self.radius
        v.norm() - self.radius
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

struct Torus {
    r1: f32,
    r2: f32,
}

impl Sdf for Torus {
    fn sdf(&self, v: &Vector3<f32>) -> f32 {
        let q: Vector2<f32> = Vector2::new(Vector2::new(v[0], v[2]).norm() - self.r1, v[1]);
        q.norm() - self.r2
    }
}


struct Scene {
    camera: Camera,
    sdf: fn(&Vector3<f32>) -> f32,
}

impl Scene {
    fn render(&self) -> image::ImageBuffer<image::Rgb<f32>, Vec<f32>> {
        let back_color: image::Rgb<f32> = image::Rgb([0., 0., 0.]);

        let mut imgbuf = image::ImageBuffer::new(self.camera.res.0 as u32, self.camera.res.1 as u32);

        let clock = std::time::Instant::now();

        for y in 0..self.camera.res.1 {

            for x in 0..self.camera.res.0 {
                let r = self.camera.get_ray(x, y);
                let y = self.camera.res.1 - y -1;  // flip vertficaly

                let dist = self.march(&r);
                let pos = r.org + *r.dir * dist;

                if dist == std::f32::INFINITY {
                    *imgbuf.get_pixel_mut(x as u32, y as u32) = back_color;
                }
                else {
                    if dist < 0. {
                        *imgbuf.get_pixel_mut(x as u32, y as u32) = image::Rgb([255., 0., 0.]);
                    }
                    else {
                        let color = [
                            pos.x.fract() * 256f32,
                            pos.y.fract() * 256f32,
                            pos.z.fract() * 256f32,

                        ];

                        *imgbuf.get_pixel_mut(x as u32, y as u32) = image::Rgb(color);
                    }
                }

                if x < 10 && y < 10 {
                    *imgbuf.get_pixel_mut(x as u32, y as u32) = image::Rgb([250., 250., 250.]);
                }

            }
        }
        let elapsed = clock.elapsed();
        println!("render time: {:?}", elapsed);

        imgbuf
    }

    fn march(&self, r: &Ray) -> f32 {
        let mut current_pos = r.org;
        let mut current_d = 0f32;

        for i in 0..MAX_STEPS {
            let d = (self.sdf)(&current_pos);

            if d < PRECISION {
                return current_d;
            }

            current_d += d;
            current_pos += *r.dir * d;
        }

        std::f32::INFINITY
    }
}

fn dist(v: &Vector3<f32>) -> f32 {
    let sphere = Sphere { center: Vector3::new(10., 20., 0.), radius: 10f32 };
    let plane = Plane { normal: Unit::new_normalize(Vector3::new(0., 1., 0.)), dist: 0. };
    let torus = Torus { r1: 4f32, r2: 1f32 };

    //plane.sdf(v).min(sphere.sdf(v))
    //sphere.sdf(&v)


    let q = Vector3::new(v.x.rem_euclid(10.), v.y, v.z.rem_euclid(10.)) - 0.5f32 * Vector3::new(10.0, 0.0, 10.);
    torus.sdf(&q)
}


fn tone_map(img: &image::ImageBuffer<image::Rgb<f32>, Vec<f32>> ) -> image::ImageBuffer<image::Rgb<u8>, Vec<u8>> {
    let clock = std::time::Instant::now();

    let n = img.width() * img.height();
    
    let mean: f32 = img.pixels().map(|p| p[0].max(p[1]).max(p[2])).sum::<f32>() / n as f32;
    let sqr_mean: f32 = img.pixels().map(|p| p[0].max(p[1]).max(p[2]) * p[0].max(p[1]).max(p[2])).sum::<f32>() / n as f32;

    let variance = sqr_mean - mean * mean;
    let max_exposure = mean + variance.sqrt();


    let mut result: image::ImageBuffer<image::Rgb<u8>, Vec<u8>> = image::ImageBuffer::new(img.width(), img.width());

    for (p_out, p_in) in result.pixels_mut().zip(img.pixels()) {
        for i in 0..3 {
            let ref mut sub_out = p_out[i];
            let sub_in = p_in[i];
            *sub_out = (sub_in.max(256f32)) as u8;
            //*sub_out = (((sub_in as f32 / max_exposure) * 256f32).max(256f32)) as u8;
        }
    }

    let elapsed = clock.elapsed();
    println!("tone mapping time: {:?}", elapsed);

    result
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

    let scene = Scene {
        camera: c,
        sdf: dist,
    };

    let raw_image = scene.render();
    let imgbuf = tone_map(&raw_image);
    image::ImageRgb8(imgbuf).save("out.png").unwrap();

}
