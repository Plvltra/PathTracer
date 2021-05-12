// a Path Tracer by Zhengxun Zhang, 2021
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <time.h>
#include <iostream>
#include <string>
using namespace std;

double M_PI = 3.14159265358979;
double M_1_PI = 1.0 / M_PI;
double erand48(unsigned short xsubi[3])
{
	return (double)rand() / (double)RAND_MAX;
}

double random(double upper = 1.0)
{
	return (double)rand() / (double)RAND_MAX * upper;
}

// Vec structure acts as points, colors, vectors
struct Vec
{
	double x, y, z;   // position, also color (r,g,b)

	Vec(double x_ = 0, double y_ = 0, double z_ = 0)
	{
		x = x_; y = y_; z = z_;
	}
	Vec operator+(const Vec& b) const
	{
		return Vec(x + b.x, y + b.y, z + b.z);
	}
	Vec operator-(const Vec& b) const
	{
		return Vec(x - b.x, y - b.y, z - b.z);
	}
	Vec operator*(double b) const
	{
		return Vec(x * b, y * b, z * b);
	}
	Vec operator/(double b) const
	{
		return Vec(x / b, y / b, z / b);
	}
	Vec mult(const Vec& b) const
	{
		return Vec(x * b.x, y * b.y, z * b.z);
	}
	Vec& norm()
	{
		return *this = *this * (1 / sqrt(x * x + y * y + z * z));
	}
	double dot(const Vec& b) const
	{
		return x * b.x + y * b.y + z * b.z;
	}
	// cross
	Vec operator%(Vec& b)
	{
		return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
	}
	double length()
	{
		return sqrt(x * x + y * y + z * z);
	}
	string str()
	{
		return to_string(x) + " " + to_string(y) + " " + to_string(z);
	}
};

// Ray structure
struct Ray
{
	Vec o, d;
	Ray(Vec o_, Vec d_) : o(o_), d(d_)
	{
	}
};

// enum of materials types used in radiance function
enum Refl_t
{
	DIFF, SPEC, REFR
};

// smallpt only supports spheres
struct Sphere
{
	double rad;     // radius
	Vec p, e, c;    // position, emission, color
	Refl_t refl;    // reflection type (DIFFuse, SPECular, REFRactive)

	// constructor: radius, center, emission, color, reflection type
	Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_) :
		rad(rad_), p(p_), e(e_), c(c_), refl(refl_)
	{
	}

	// 0: miss hit
	double intersect(const Ray& r) const
	{
		// solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
		Vec op = p - r.o;                       // p is sphere center (C)
		double t, eps = 1e-4;                   // eps is a small fudge factor
		double b = op.dot(r.d);                 // 1/2 b from quadratic eq. setup
		double det = b * b - op.dot(op) + rad * rad;    // (b^2-4ac)/4: a=1 because ray normalised
		if (det < 0) return 0;                  // ray misses sphere
		else det = sqrt(det);
		return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0); // return smaller positive t
	}
};

// hard coded scene description
// the scene description consists of a bunch of spheres
// scene: radius, position, color, material

Sphere spheres[] = {
	Sphere(1e5, Vec(1e5 + 1,40.8,81.6),	 Vec(),				Vec(.75,.25,.25),	DIFF), //Left
	Sphere(1e5, Vec(-1e5 + 99,40.8,81.6),Vec(),				Vec(.25,.25,.75),	DIFF), //Rght
	Sphere(1e5, Vec(50,40.8, 1e5),		 Vec(),				Vec(.75,.75,.75),	DIFF), //Back
	Sphere(1e5, Vec(50,40.8,-1e5 + 170), Vec(),				Vec(),				DIFF), //Frnt
	Sphere(1e5, Vec(50, 1e5, 81.6),		 Vec(),				Vec(.75,.75,.75),	DIFF), //Botm
	Sphere(1e5, Vec(50,-1e5 + 81.6,81.6),Vec(),				Vec(.75,.75,.75),	DIFF), //Top
	Sphere(16.5,Vec(27,16.5,47),		 Vec(),				Vec(1,1,1) * .999,	SPEC), //Mirr
	Sphere(16.5,Vec(73,16.5,78),		 Vec(),				Vec(1,1,1) * .999,	REFR), //Glas
	Sphere(1.5, Vec(50,81.6 - 16.5,81.6),Vec(4,4,4) * 100,	Vec(),				DIFF)  //Lite
};
int numSpheres = sizeof(spheres) / sizeof(Sphere);
int lightID = numSpheres - 1;

// clamp x in [0, 1]
inline double clamp(double x)
{
	return x < 0 ? 0 : x>1 ? 1 : x;
}

// converts floats to integers to be saved in ppm file
inline int toInt(double x)
{
	return int(pow(clamp(x), 1 / 2.2) * 255 + .5);
}

// intersects ray with scene
inline bool intersect(const Ray& r, double& t, int& id)
{
	double n = sizeof(spheres) / sizeof(Sphere);
	double d;
	double inf = t = 1e20;

	for (int i = int(n); i--;) {
		if ((d = spheres[i].intersect(r)) && d < t) {
			t = d;
			id = i;
		}
	}
	return t < inf;
}

#define _ZZX_
#ifdef _ZZX_
int sample_count = 2; 
const double P_RR = 0.5;
// TODO: fix pdf, fr
Vec calcPixelColor(Ray ray, int levels)
{
	double t;	// hit distance
	int id;		// hit id
	bool hit = intersect(ray, t, id);
	if (!hit)
		return Vec(); // black background

	Sphere obj = spheres[id];
	Vec p = ray.o + ray.d * t;		// hit point
	Vec n = (p - obj.p).norm();		// sphere's surface normal
	Vec nl = ray.d.dot(n) < 0 ? n : n * -1; // real normal(normal will points to inner if ray is in glass)
	
	// direct light
	if (obj.refl == DIFF) 
	{
		Vec L_dir = 0.0;
		double radius = spheres[lightID].rad;
		Vec sample = spheres[lightID].p + Vec(random(radius), random(radius), random(radius)); // random point in light
		Ray light_ray = Ray(p, (sample - p).norm());
		double t1;
		int id1;
		if (intersect(light_ray, t1, id1) && lightID == id1)
		{
			Vec x_p = light_ray.o + light_ray.d * t1;
			Vec n_p = (x_p - spheres[lightID].p).norm();// hit normal
			double cos_theta = fabs(light_ray.d.dot(n));// TODO:
			double cos_theta_p = -light_ray.d.dot(n_p);
			double dist = (x_p - p).length();
			double pdf_light = 1 / (1.5 * 1.5 * 1.5);
			L_dir = spheres[lightID].e.mult(obj.c) * cos_theta * cos_theta_p / (dist * dist) / pdf_light;
		}

		if (random() > P_RR) {
			return L_dir;
		}

		// indirect light
		double r1 = 2 * M_PI * rand();  // angle around
		double r2 = random(), r2s = sqrt(r2);  // distance from center
		Vec w = n;  // w = normal
		Vec u = (ray.d - n * ray.d.dot(n)).norm(); // u is perpendicular to w
		Vec v = w % u;  // v is perpendicular to u and w
		Vec refl_dir = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm(); // reflect dir
		Vec L_indir = 0.0;
		double t2;
		int id2;
		Ray refl_ray = Ray(p, refl_dir);
		if (intersect(refl_ray, t2, id2))
		{
			Vec fr = obj.c;
			double cos_theta = n.dot(refl_dir), pdf_hemi = 1 / (2 * M_PI);
			L_indir = calcPixelColor(refl_ray, levels + 1).mult(fr) * cos_theta / pdf_hemi / P_RR;
		}
		return obj.e + L_dir + L_indir;
	}	
	else if (obj.refl == SPEC)
	{
		Vec refl = (ray.d - n * 2 * n.dot(ray.d)).norm(); // reflect dir
		Vec fr = obj.c;
		double cos_theta = n.dot(refl);
		return obj.e + calcPixelColor(Ray(p, refl), levels + 1).mult(fr) * cos_theta;
		//TODO: recur end
	}
	else if (obj.refl == REFR)
	{
		if (levels > 3)
			return Vec();

		Vec fr = obj.c;
		Vec refl = ray.d - nl * 2 * nl.dot(ray.d); // reflect dir
		bool into = n.dot(nl) > 0;
		double nc = 1, nt = 1.5; // refraction values of air and glass
		double nnt = into ? nc / nt : nt / nc; // incident dielectric / refraction dielectric
		double ddn = ray.d.dot(nl);
		double cos2t = 1 - nnt * nnt * (1 - ddn * ddn);
		// TOTAL REFLECTION OCCURS
		if (cos2t < 0){
			return obj.e + calcPixelColor(Ray(p, refl), levels + 1).mult(fr);
		}

		Vec refr = (ray.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm(); // refract dir
		// calculate fresnel term
		double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : refr.dot(n));
		//double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
		double Re = R0 + (1 - R0) * c * c * c * c * c; // Re: reflectance
		double Tr = 1 - Re;
		return obj.e + calcPixelColor(Ray(p, refl), levels + 1).mult(fr) * Re
			+ calcPixelColor(Ray(p, refr), levels + 1).mult(fr) * Tr;
	}
}

// main function, loops over image pixels, creates image,
// and saves it to a ppm file
int main(int argc, char* argv[])
{
	srand((unsigned)time(0));
	int w = 512, h = 384;  // image size
	int samps = argc == 2 ? atoi(argv[1]) / 4 : 1; // # samples (default of 1)
	Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // cam pos, dir
	Vec cx = Vec(w * .5135 / h);  // x direction increment (uses implicit 0 for y, z)
	Vec cy = (cx % cam.d).norm() * .5135;  // y direction increment (% : cross product)
	Vec r;  // used for colors of samples
	Vec* c = new Vec[w * h];  // the image

	for (int y = 0; y < h; y++) 
	{
		fprintf(stderr, "\rRendering %5.2f%%", 100.0 * y / h);
		for(int x = 0; x < w; x++)
		{
			int idx = x + y * w;
			double gap = 1.0 / (sample_count + 1);
			for (int sx = 0; sx < sample_count; sx++) // add a little bias to sample point
				for (int sy = 0; sy < sample_count; sy++) {
					double px = x - 0.5 + sx * gap;
					double py = h - 1 - (y - 0.5 + sy * gap);
					// right-handed-coordinate
					Vec rd = cx * (px / w - 0.5) + cy * (py / h - 0.5) + cam.d;
					rd.norm();
					Vec ro = cam.o + rd * 140;
					c[idx] = c[idx] + calcPixelColor(Ray(ro, rd), 0);
				}
			c[idx] = c[idx] / (sample_count * sample_count);
		}
	}

	// write out the file to a ppm
	FILE* f;
	fopen_s(&f, "image.ppm", "w");         // Write image to PPM file.
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
	for (int i = 0; i < w * h; i++)
		fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
	return 0;
}

#else
// compute the radiance estimate along ray r
Vec radiance(const Ray& r, int depth, unsigned short* Xi, int E = 1)
{
	double t;                               // distance to intersection
	int id = 0;                               // id of intersected object
	if (!intersect(r, t, id)) return Vec(); // if miss, return black
	const Sphere& obj = spheres[id];        // the hit object

	if (depth > 10) return Vec();

	Vec x = r.o + r.d * t;  // ray intersection point
	Vec n = (x - obj.p).norm();  // sphere normal
	Vec nl = n.dot(r.d) < 0 ? n : n * -1;  // properly oriented surface normal
	Vec f = obj.c;  // object color (BRDF modulator)
// use maximum reflectivity amount of Russian roulette
	double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; // max refl
	if (++depth > 5 || !p) {
		if (erand48(Xi) < p)
			f = f * (1 / p);
		else
			return obj.e * E; //R.R.
	}

// ideal difuse reflection
	if (obj.refl == DIFF) {        // ideal DIFFUSE reflection
		double r1 = 2 * M_PI * erand48(Xi);  // angle around
		double r2 = erand48(Xi), r2s = sqrt(r2);  // distance from center
		Vec w = nl;  // w = normal
		Vec u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm(); // u is perpendicular to w
		Vec v = w % u;  // v is perpendicular to u and w
		Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();  // d is random reflection ray
	// loop over any lights (explicit lighting)
		Vec e;
		for (int i = 0; i < numSpheres; i++) {
			const Sphere& s = spheres[i];
			if (s.e.x <= 0 && s.e.y <= 0 && s.e.z <= 0) continue;  // skip non-lights
	// create random direction towards sphere using method from Realistic ray Tracing
			Vec sw = s.p - x, su = ((fabs(sw.x) > .1 ? Vec(0, 1) : Vec(1)) % sw).norm(), sv = sw % su;
			double cos_a_max = sqrt(1 - s.rad * s.rad / (x - s.p).dot(x - s.p));
			double eps1 = erand48(Xi), eps2 = erand48(Xi);
			double cos_a = 1 - eps1 + eps1 * cos_a_max;
			double sin_a = sqrt(1 - cos_a * cos_a);
			double phi = 2 * M_PI * eps2;
			Vec l = su * cos(phi) * sin_a + sv * sin(phi) * sin_a + sw * cos_a;
			l.norm();

			// shoot shadow ray
			if (intersect(Ray(x, l), t, id) && id == i) {  // shadow ray
				double omega = 2 * M_PI * (1 - cos_a_max);
				e = e + f.mult(s.e * l.dot(nl) * omega) * M_1_PI;  // 1/pi for brdf
			}
		}
		return obj.e * E + e + f.mult(radiance(Ray(x, d), depth, Xi, 0));
		// ideal specular reflection
	} else if (obj.refl == SPEC) {            // Ideal SPECULAR reflection
		return obj.e + f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth, Xi));
		// otherwise we have a dielectric (glass) surface
	} else {
		Ray reflRay(x, r.d - n * 2 * n.dot(r.d));     // Ideal dielectric REFRACTION
		bool into = n.dot(nl) > 0;                // Ray from outside going in?
		double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t;
		// if total internal reflection, reflect
		if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)    // Total internal reflection
			return obj.e + f.mult(radiance(reflRay, depth, Xi));
		// otherwise, choose reflection or refraction
		Vec tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
		double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : tdir.dot(n));
		double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
		return obj.e + f.mult(depth > 2 ? (erand48(Xi) < P ?   // Russian roulette
			radiance(reflRay, depth, Xi) * RP : radiance(Ray(x, tdir), depth, Xi) * TP) :
			radiance(reflRay, depth, Xi) * Re + radiance(Ray(x, tdir), depth, Xi) * Tr);
	}
}

// main function, loops over image pixels, creates image,
// and saves it to a ppm file
int main(int argc, char* argv[])
{
	int w = 512, h = 384;  // image size
	int samps = argc == 2 ? atoi(argv[1]) / 4 : 1; // # samples (default of 1)
	Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // cam pos, dir
	Vec cx = Vec(w * .5135 / h);  // x direction increment (uses implicit 0 for y, z)
	Vec cy = (cx % cam.d).norm() * .5135;  // y direction increment (% : cross product)
	Vec r;  // used for colors of samples
	Vec* c = new Vec[w * h];  // the image
#pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP
// loop over all image pixels
	for (int y = 0; y < h; y++) {                       // Loop over image rows
		fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100. * y / (h - 1));  // print progress
		unsigned short Xi[3] = { 0,0,y * y * y };
		for (unsigned short x = 0; x < w; x++)          // Loop cols
// for each pixel do 2x2 subsamples and samps samples per subsample
for (int sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++)     // 2x2 subpixel rows
	for (int sx = 0; sx < 2; sx++, r = Vec()) {        // 2x2 subpixel cols
		for (int s = 0; s < samps; s++) {
			// I believe yhis produce a tent filter
			double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
			double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
			Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
				cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
			Ray ray = Ray(cam.o + d * 140, d.norm());
			r = r + radiance(ray, 0, Xi) * (1. / samps);
		} // Camera rays are pushed ^^^^^ forward to start in interior
		c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
	}
	}
	// write out the file to a ppm
	FILE* f;
	fopen_s(&f, "image.ppm", "w");         // Write image to PPM file.
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
	for (int i = 0; i < w * h; i++)
		fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
	return 0;
}
#endif