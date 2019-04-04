#include "bsdf.h"

#include <iostream>
#include <algorithm>
#include <utility>

using std::min;
using std::max;
using std::swap;

namespace CGL {

void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {
  Vector3D z = Vector3D(n.x, n.y, n.z);
  Vector3D h = z;
  if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z)) h.x = 1.0;
  else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z)) h.y = 1.0;
  else h.z = 1.0;

  z.normalize();
  Vector3D y = cross(h, z);
  y.normalize();
  Vector3D x = cross(z, y);
  x.normalize();

  o2w[0] = x;
  o2w[1] = y;
  o2w[2] = z;
}

// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
	// TODO: 1.2
	// Using BSDF::reflect(), implement sample_f for a mirror surface
	*pdf = 1.0;
	reflect(wo, wi);
	return reflectance / abs_cos_theta(*wi);
}

// Microfacet BSDF //

double MicrofacetBSDF::G(const Vector3D& wo, const Vector3D& wi) {
    return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D& h) {
	// TODO: 2.2
	// Compute Beckmann normal distribution function (NDF) here.
	// You will need the roughness alpha.
	//std::cout << h.norm() << std::endl;
	//std::cout << alpha << std::endl;
	//double thetah = acos(h.z);
	//double thetah = getTheta(h);
	
	//if (thetah < 0) {
	//	std::cout << thetah << std::endl;
	//}
	//std::cout << pow(sin_theta(h) / cos_theta(h), 2) / pow(alpha, 2) << std::endl;
	double Dnum = exp(-pow(sin_theta(h)/cos_theta(h), 2) / pow(alpha, 2));
	double Dden = PI * pow(alpha, 2) * pow(h.z, 4);
	//double Dden = PI * pow(alpha, 2) * pow(abs_cos_theta(h), 4);
	//std::cout << Dnum / Dden << std::endl;
	return Dnum / Dden;
}

Spectrum MicrofacetBSDF::F(const Vector3D& wi) {
	// TODO: 2.3
	// Compute Fresnel term for reflection on dielectric-conductor interface.
	// You will need both eta and etaK, both of which are Spectrum.
	//std::cout << wi.norm() << std::endl;
	double Rslastterm = pow(wi.z, 2);
	Spectrum Rslast = Spectrum(Rslastterm, Rslastterm, Rslastterm);
	Spectrum Rsnum = (eta*eta) + (k*k) - (eta*2.0*wi.z) + Rslast;
	Spectrum Rsden = (eta*eta) + (k*k) + (eta*2.0*wi.z) + Rslast;

	double Rplastterm = 1.0;
	Spectrum Rplast = Spectrum(Rplastterm, Rplastterm, Rplastterm);
	Spectrum Rpnum = ((eta*eta) + (k*k))*pow(wi.z, 2) - (eta*2.0*wi.z) + Rplast;
	Spectrum Rpden = ((eta*eta) + (k*k))*pow(wi.z, 2) + (eta*2.0*wi.z) + Rplast;

	Spectrum Rs = Rsnum / Rsden;
	Spectrum Rp = Rpnum / Rpden;

	return (Rs + Rp) / 2.0;
}

Spectrum MicrofacetBSDF::f(const Vector3D& wo, const Vector3D& wi) {
	// TODO: 2.1
	// Implement microfacet model here
	//std::cout << wi.norm() << std::endl;
	if (!(wo.z > 0 && wi.z > 0)) return Spectrum();
	Vector3D wowi = wo + wi;
	Vector3D h = wowi / wowi.norm();
	//std::cout << h.norm() << std::endl;
	//return Spectrum();
	return F(wi) * G(wo, wi) * D(h) / (4.0 * wo.z * wi.z);
}

Spectrum MicrofacetBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
	// TODO: 2.4
	// *Importance* sample Beckmann normal distribution function (NDF) here.
	// Note: You should fill in the sampled direction *wi and the corresponding *pdf,
	//       and return the sampled BRDF value.
	//std::cout << alpha << std::endl;
	//*wi = cosineHemisphereSampler.get_sample(pdf); //placeholder
	//return MicrofacetBSDF::f(wo, *wi);
	//std::cout << wo.norm() << std::endl;
	Vector2D r = sampler.get_sample(); // r1, r2 in [0, 1)

	float thetah = atan(sqrt(-pow(alpha, 2) * log(1.0 - r[0])));
	float phih = 2 * PI * r[1];

	float ptheta = 2.0 * sin(thetah) / (pow(alpha, 2) * pow(cos(thetah), 3)) * exp(-pow(tan(thetah), 2) / pow(alpha, 2));
	float pphi = 1.0 / (2.0 * PI);

	Vector3D h = Vector3D(sin(thetah)*cos(phih), sin(thetah)*sin(phih), cos(thetah));

	*wi = -wo + 2.0*dot(wo, h)*h;
	if (!wi->z > 0) {
		*pdf = 0;
		return Spectrum();
	}
	//std::cout << wi->norm() << std::endl;

	float pwh = ptheta * pphi / sin(thetah);
	*pdf = pwh / (4.0 * dot(*wi, h));

	return MicrofacetBSDF::f(wo, *wi);
}

// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  return Spectrum();
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
	// TODO: 1.4
	// Compute Fresnel coefficient and either reflect or refract based on it.
	bool tir = !refract(wo, wi, ior);
	if (tir) { // total internal reflection
		reflect(wo, wi);
		*pdf = 1.0;
		return reflectance / abs_cos_theta(*wi);
	}
	else {
		float R0 = pow((1.0 - ior)/(1.0 + ior), 2);
		float R = R0 + (1.0 - R0)*pow(1.0 - abs(wo.z), 5);
		if (coin_flip(R)) {
			reflect(wo, wi);
			*pdf = R;
			return R * reflectance / abs_cos_theta(*wi);
		}
		else {
			float eta = (wo.z < 0) ? ior : 1.0 / ior;
			*pdf = 1.0 - R;
			return (1.0 - R) * transmittance / abs_cos_theta(*wi) / pow(eta, 2);
		}
	}
}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {
	// TODO: 1.1
	// Implement reflection of wo about normal (0,0,1) and store result in wi.
	Vector3D n(0, 0, 1);
	*wi = -wo + 2.0 * wo.z * n;
}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) { // ior is ratio of old ior to new ior
	// TODO: 1.3
	// Use Snell's Law to refract wo surface and store result ray in wi.
	// Return false if refraction does not occur due to total internal reflection
	// and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
	// ray entering the surface through vacuum.
	float sign_woz = wo.z >= 0 ? 1.0 : -1.0;
	
	float eta = (wo.z < 0) ? ior : 1.0 / ior;

	float costhetap = 1.0 - pow(eta, 2)*(1.0 - pow(wo.z, 2));
	if (costhetap < 0) return false; // total internal reflection

	*wi = Vector3D(-wo.x * eta, -wo.y * eta, -sign_woz * sqrt(costhetap));
	return true;
}

// Emission BSDF //

Spectrum EmissionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1.0 / PI;
  *wi  = sampler.get_sample(pdf);
  return Spectrum();
}

} // namespace CGL
