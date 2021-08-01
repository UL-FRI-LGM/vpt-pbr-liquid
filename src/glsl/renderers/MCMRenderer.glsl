// #package glsl/shaders

// #include ../mixins/Photon.glsl
// #include ../mixins/rand.glsl
// #include ../mixins/unprojectRand.glsl
// #include ../mixins/intersectCube.glsl

// #section MCMGenerate/vertex

void main() {}

// #section MCMGenerate/fragment

void main() {}

// #section MCMIntegrate/vertex

#version 300 es

layout (location = 0) in vec2 aPosition;

out vec2 vPosition;

void main() {
    vPosition = aPosition;
    gl_Position = vec4(aPosition, 0.0, 1.0);
}

// #section MCMIntegrate/fragment

#version 300 es
precision mediump float;

#define M_INVPI 0.31830988618
#define M_2PI 6.28318530718
#define EPS 1e-5
#define SQRT3 1.73205080757

#define MSIZE 255

#define AVOGADRO_CONSTANT 6.02214179 * pow(10.0f, 23.0f)
#define MOLAR_MASS_WATER 18.01528 * pow(10.0f, -3.0f)
#define M_PI 3.1415926535
#define POLARIZABILITY_WATER_532nm 1.4864 * pow(10.0f, 30.0f)
#define FLOOR_DENSITY 3000.0f
#define CUBE_DENSTIY 2500.0f
#define FLOOR_REFRACTIVE_INDEX 50000.0f
#define CUBE_REFRACTIVE_INDEX 30000.0f

@Photon

uniform mediump sampler2D uPosition;
uniform mediump sampler2D uDirection;
uniform mediump sampler2D uTransmittance;
uniform mediump sampler2D uRadiance;

uniform mediump sampler2D uEnvironment;

uniform mediump sampler2D uTransferFunction0;
uniform mediump sampler3D uVolume0;

uniform mediump sampler2D uTransferFunction1;
uniform mediump sampler3D uVolume1;

uniform mediump sampler2D uTransferFunction2;
uniform mediump sampler3D uVolume2;

uniform mediump sampler2D uTransferFunction3;
uniform mediump sampler3D uVolume3;

uniform int uNumberOfChannels;
uniform vec4 uChannelContributions;

uniform mat4 uMvpInverseMatrix;
uniform vec2 uInverseResolution;
uniform float uRandSeed;
uniform float uBlur;

uniform float uAbsorptionCoefficient;
uniform float uScatteringCoefficient;
uniform float uScatteringBias;
uniform float uMajorant;
uniform uint uMaxBounces;
uniform uint uSteps;
uniform bool uMaxContribution;
uniform bool uOrigData;
uniform float uOrigVsSeg;
uniform bool uBilateral;
uniform bool uBilateralGradient;
// Slow implementation. For improvement, put this variables int define
uniform float uBilateralSIGMA; // 10.0
uniform float uBilateralBSIGMA; // 0.1
uniform int uBilateralMSIZE; // 15

uniform mat4 uEnvironmentRotationMatrix;
uniform bool uEnvironmentTextureOverride;
uniform vec3 uEnvironmentColor;
uniform float uEnvironmentContribution;

uniform vec3 uMinCutPlaneValues;
uniform vec3 uMaxCutPlaneValues;
uniform float uViewCutDistance;

uniform float uMinDensity;
uniform float uMaxDensity;

uniform float uSize;

in vec2 vPosition;

layout (location = 0) out vec4 oPosition;
layout (location = 1) out vec4 oDirection;
layout (location = 2) out vec4 oTransmittance;
layout (location = 3) out vec4 oRadiance;

@rand
@unprojectRand
@intersectCube

float normpdf(in float x, in float sigma)
{
  return 0.39894*exp(-0.5*x*x/(sigma*sigma))/sigma;
}

float normpdf2(in vec2 v, in float sigma)
{
  return 0.39894*exp(-0.5*dot(v,v)/(sigma*sigma))/sigma;
}

float normpdf3(in vec3 v, in float sigma)
{
  return 0.39894*exp(-0.5*dot(v,v)/(sigma*sigma))/sigma;
}

float bilateralFiltering3D(vec3 position, float volumeSample) {
    //declare stuff
    int kSize = (uBilateralMSIZE-1)/2;
    float kernel[MSIZE];
    float bfinal_colour = 0.0;

    float bZ = 0.0;

    //create the 1-D kernel
    for (int j = 0; j <= kSize; ++j) {
      kernel[kSize+j] = kernel[kSize-j] = normpdf(float(j), uBilateralSIGMA);
    }

    float cc;
    float gfactor;
    float bfactor;
    float bZnorm = 1.0/normpdf(0.0, uBilateralBSIGMA);
    //read out the texels
    for (int i=-kSize; i <= kSize; ++i) {
        for (int j=-kSize; j <= kSize; ++j) {
            for (int k=-kSize; k <= kSize; ++k) {
                // voxel values in the 2D neighborhood
                vec3 coord = position.xyz + vec3(float(i), float(j), float(k));
                cc = texture(uVolume0, coord).r;

                // compute both the gaussian smoothed and bilateral
                gfactor = kernel[kSize+k]*kernel[kSize+j]*kernel[kSize+i];
                bfactor = normpdf(cc-volumeSample, uBilateralBSIGMA)*bZnorm*gfactor;
                bZ += bfactor;

                bfinal_colour += bfactor*cc;
            }
        }
    }
    return bfinal_colour/bZ;
}

vec2 bilateralFiltering3D_2(vec3 position, vec2 volumeSample) {
    //declare stuff
    int kSize = (uBilateralMSIZE-1)/2;
    float kernel[MSIZE];
    vec2 bfinal_colour = vec2(0.0);

    float bZ = 0.0;

    //create the 1-D kernel
    for (int j = 0; j <= kSize; ++j) {
      kernel[kSize+j] = kernel[kSize-j] = normpdf(float(j), uBilateralSIGMA);
    }

    vec2 cc;
    float gfactor;
    float bfactor;
    float bZnorm = 1.0/normpdf(0.0, uBilateralBSIGMA);
    //read out the texels
    for (int i=-kSize; i <= kSize; ++i) {
        for (int j=-kSize; j <= kSize; ++j) {
            for (int k=-kSize; k <= kSize; ++k) {
                // voxel values in the 2D neighborhood
                vec3 coord = position.xyz + vec3(float(i), float(j), float(k));
                cc = texture(uVolume0, coord).rg;

                // compute both the gaussian smoothed and bilateral
                gfactor = kernel[kSize+k]*kernel[kSize+j]*kernel[kSize+i];
                bfactor = normpdf2(cc-volumeSample, uBilateralBSIGMA)*bZnorm*gfactor;
                bZ += bfactor;

                bfinal_colour += bfactor*cc;
            }
        }
    }
    return bfinal_colour/bZ;
}

void resetPhoton(inout vec2 randState, inout Photon photon) {
    vec3 from, to;
    unprojectRand(randState, vPosition, uMvpInverseMatrix, uInverseResolution, uBlur, from, to);
    photon.direction = normalize(to - from);
    photon.bounces = 0u;
    vec2 tbounds = max(limitedIntersectCube(from, photon.direction, uMinCutPlaneValues, uMaxCutPlaneValues), 0.0);
    photon.position = from + (tbounds.x + uViewCutDistance * SQRT3) * photon.direction;
    // photon.position = from + tbounds.x * photon.direction;
    photon.transmittance = vec3(1);
}

vec4 sampleEnvironmentMap(vec3 d) {
    vec4 rotatedD = uEnvironmentRotationMatrix * vec4(d, 1.0);
    vec2 texCoord = vec2(atan(rotatedD.x, -rotatedD.z), asin(-rotatedD.y) * 2.0) * M_INVPI * 0.5 + 0.5;
    // vec2 texCoord = vec2(atan(d.x, -d.z), asin(-d.y) * 2.0) * M_INVPI * 0.5 + 0.5;
    return texture(uEnvironment, texCoord);
}

vec4 sampleVolumeColor(vec3 position) {
    vec4 channelContribs = uChannelContributions / max(max(max(uChannelContributions.x, uChannelContributions.y), uChannelContributions.z), uChannelContributions.w);

    vec2 volumeSample0 = texture(uVolume0, position).rg ;
    vec4 transferSample0 = texture(uTransferFunction0, volumeSample0) * channelContribs.x;
    vec4 transferSample1 = vec4(0);
    vec4 transferSample2 = vec4(0);
    vec4 transferSample3 = vec4(0);
    float maxVolValue = volumeSample0.r;
    vec4 maxVolSample = transferSample0;

    if (uBilateral && !uBilateralGradient) {
        volumeSample0.r = bilateralFiltering3D(position, volumeSample0.r);
    } else if (uBilateral && uBilateralGradient) {
        volumeSample0.rg = bilateralFiltering3D_2(position, volumeSample0);
    }

    if (uNumberOfChannels > 1) {
        vec2 volumeSample1 = texture(uVolume1, position).rg ;
        transferSample1 = texture(uTransferFunction1, volumeSample1) * channelContribs.y;

        if (uBilateral && !uBilateralGradient) {
            volumeSample1.r = bilateralFiltering3D(position, volumeSample1.r);
        } else if (uBilateral && uBilateralGradient) {
            volumeSample1.rg = bilateralFiltering3D_2(position, volumeSample1);
        }

        maxVolValue = volumeSample1.r;
        maxVolSample = transferSample1;
        if (maxVolValue < volumeSample1.r) {
             maxVolValue = volumeSample1.r;
             maxVolSample = transferSample1;
        }
    }
    if (uNumberOfChannels > 2) {
        vec2 volumeSample2 = texture(uVolume2, position).rg ;
        transferSample2 = texture(uTransferFunction2, volumeSample2) * channelContribs.z;

        if (uBilateral && !uBilateralGradient) {
            volumeSample2.r = bilateralFiltering3D(position, volumeSample2.r);
        } else if (uBilateral && uBilateralGradient) {
            volumeSample2.rg = bilateralFiltering3D_2(position, volumeSample2);
        }

        if (maxVolValue < volumeSample2.r) {
             maxVolValue = volumeSample2.r;
             maxVolSample = transferSample2;
        }
    }
    if (uNumberOfChannels > 3) {
        vec2 volumeSample3 = texture(uVolume3, position).rg ;
        transferSample3 = texture(uTransferFunction3, volumeSample3) * channelContribs.w;

        if (uBilateral && !uBilateralGradient) {
            volumeSample3.r = bilateralFiltering3D(position, volumeSample3.r);
        } else if (uBilateral && uBilateralGradient) {
            volumeSample3.rg = bilateralFiltering3D_2(position, volumeSample3);
        }

        if (maxVolValue < volumeSample3.r) {
             maxVolValue = volumeSample3.r;
             maxVolSample = transferSample3;
        }
    }

    if (!uMaxContribution) {
        float sumA = (transferSample0.a + transferSample1.a + transferSample2.a + transferSample3.a);
        vec3 sumC = (transferSample0.rgb * transferSample0.a + transferSample1.rgb * transferSample1.a + transferSample2.rgb * transferSample2.a + transferSample3.rgb * transferSample3.a) / sumA;
        return vec4(sumC, sumA / float(uNumberOfChannels));
    } else {
        if (uOrigData) {
            float alpha = uOrigVsSeg * transferSample0.a + (1.0 - uOrigVsSeg) * maxVolSample.a;
            vec3 color = (uOrigVsSeg * transferSample0.rgb * transferSample0.a + (1.0 - uOrigVsSeg) * maxVolSample.rgb * maxVolSample.a);
            return vec4(color, alpha);
        } else {
            return maxVolSample;
        }
    }
}

// calculate actual density value from .raw file
float calculateDensityFromRatio(float ratio) {
    if (ratio == 255.0f) {
        return FLOOR_DENSITY;
    } else if (ratio == 254.0f) {
        return CUBE_DENSTIY;
    }
    float interval = uMaxDensity - uMinDensity;
    float percentage = ratio / 255.0f;
    return (percentage * interval) + uMinDensity;
}

float calculateRefractiveIndexFromDensity(float d) {
    if (d == FLOOR_DENSITY) {
        // if density is that of floor, return very big refractive index to assure bouncing
        return FLOOR_REFRACTIVE_INDEX;
    } else if (d == CUBE_DENSTIY) {
        return CUBE_REFRACTIVE_INDEX;
    }
    float x = (4.0f * M_PI * d * AVOGADRO_CONSTANT * POLARIZABILITY_WATER_532nm / (3.0f * MOLAR_MASS_WATER)) + 1.0f;
    return pow(x, 0.66f);
}

vec3 getVectorOfCertainLength(vec3 vector, float desiredLength) {
    float vectorLength = length(vector);
    return desiredLength * vector / vectorLength;
}

float samplePhotonAndCalculateRefractiveIndex(Photon photon) {
    // samplenja fotona v smeri in izračun lomnega količnika tam
    // premik fotona v trenutni smeri za natanko velikosti voksla
    vec3 photonDirectionNextVoxel = getVectorOfCertainLength(photon.direction, 1.0f / uSize);
    vec3 newPhotonPosition = photon.position + photonDirectionNextVoxel;

    // izračun gostote in lomenga kolicnika v tem vokslu
    float mappedDensity = texture(uVolume0, newPhotonPosition).r;
    float density = calculateDensityFromRatio(mappedDensity);

    float refractiveIndex = calculateRefractiveIndexFromDensity(density);
    return refractiveIndex;
}

float calculateRefractiveAngle(float alpha, float n1, float n2) {
    // izracun lomnega kota z lomnim zakonom
    float stevec = sin(alpha) * n1;
    return asin(stevec / n2);
}

// izracun novega smernega vektorja po lomu
// todo: mozno, da bo ze zacetni cos < 0, potem resetiraj glej https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
vec3 calculateNewDirectionVector(vec3 incomingVector, vec3 normal, float riRatio, float alpha, float beta) {
    vec3 ratioIncoming = riRatio * incomingVector;
    float ratioCosine = riRatio * cos(alpha) - cos(beta);
    vec3 ratioCosineNormal = ratioCosine * normal;
    return ratioIncoming + ratioCosineNormal;
}

vec3 refractPhoton(Photon photon, vec3 gradientVector, float n1, float n2) {
    // handlanje negativnega cos, da je normala obratno obrnjena
    if (dot(gradientVector * -1.0f, photon.direction) < 0.0f)
        gradientVector = gradientVector * -1.0f;
    float alpha = acos(dot(normalize(photon.direction), normalize(gradientVector)));
    float beta = calculateRefractiveAngle(alpha, n1, n2);
    return calculateNewDirectionVector(photon.direction, gradientVector, n1 / n2, alpha, beta);
}

vec3 bouncePhoton(Photon photon, vec3 gradientVector) {
    // bounce photon
    return reflect(photon.direction, normalize(gradientVector));
}

vec3 determineNewPhotonDirection(Photon photon) {
    // calculation of current refractive index
    float currentMappedDensity = texture(uVolume0, photon.position).r;
    vec3 gradientVector = texture(uVolume0, photon.position).gba;
    float currentDensity = calculateDensityFromRatio(currentMappedDensity);
    float currentRI = calculateRefractiveIndexFromDensity(currentDensity);

    // caluclation of next refractive index
    float nextRI = samplePhotonAndCalculateRefractiveIndex(photon);

    // determine if perform bounce or refract
        // TODO: determine right mark for refraction or bounce
    if (abs(currentRI - nextRI) < 1000.0f) {
        // refract
        return refractPhoton(photon, gradientVector, currentRI, nextRI);
    } else {
        // bounce (or return the same direction)
        //return bouncePhoton(photon, gradientVector);
        return photon.direction;
    }

}

vec3 randomDirection(vec2 U) {
    float phi = U.x * M_2PI;
    float z = U.y * 2.0 - 1.0;
    float k = sqrt(1.0 - z * z);
    return vec3(k * cos(phi), k * sin(phi), z);
}

float sampleHenyeyGreensteinAngleCosine(float g, float U) {
    float g2 = g * g;
    float c = (1.0 - g2) / (1.0 - g + 2.0 * g * U);
    return (1.0 + g2 - c * c) / (2.0 * g);
}

vec3 sampleHenyeyGreenstein(float g, vec2 U, vec3 direction) {
    // generate random direction and adjust it so that the angle is HG-sampled
    vec3 u = randomDirection(U);
    if (abs(g) < EPS) {
        return u;
    }
    float hgcos = sampleHenyeyGreensteinAngleCosine(g, fract(sin(U.x * 12345.6789) + 0.816723));
    float lambda = hgcos - dot(direction, u);
    return normalize(u + lambda * direction);
}

void main() {
    Photon photon;
    vec2 mappedPosition = vPosition * 0.5 + 0.5;
    photon.position = texture(uPosition, mappedPosition).xyz;
    vec4 directionAndBounces = texture(uDirection, mappedPosition);
    float density = texture(uVolume0, photon.position).w;
    photon.density = density;
    photon.previousDensity = density;
    photon.direction = directionAndBounces.xyz;
    photon.bounces = uint(directionAndBounces.w + 0.5);
    photon.transmittance = texture(uTransmittance, mappedPosition).rgb;
    vec4 radianceAndSamples = texture(uRadiance, mappedPosition);
    photon.radiance = radianceAndSamples.rgb;
    photon.samples = uint(radianceAndSamples.w + 0.5);

    vec2 r = rand(vPosition * uRandSeed);
    for (uint i = 0u; i < uSteps; i++) {
        r = rand(r);
        float t = -log(r.x) / uMajorant;
        photon.position += t * photon.direction;

        vec4 volumeSample = sampleVolumeColor(photon.position);
        float muAbsorption = volumeSample.a * uAbsorptionCoefficient;
        float muScattering = volumeSample.a * uScatteringCoefficient;
        float muNull = uMajorant - muAbsorption - muScattering;
        float muMajorant = muAbsorption + muScattering + abs(muNull);
        float PNull = abs(muNull) / muMajorant;
        float PAbsorption = muAbsorption / muMajorant;
        float PScattering = muScattering / muMajorant;

        photon.density = length(texture(uVolume0, photon.position).rgb);

        // if (any(greaterThan(photon.position, vec3(1))) || any(lessThan(photon.position, vec3(0)))) {
        if (any(greaterThan(photon.position, uMaxCutPlaneValues)) || any(lessThan(photon.position, uMinCutPlaneValues))) {
            // out of bounds
            vec4 envSample = sampleEnvironmentMap(photon.direction);
            if (uEnvironmentTextureOverride && photon.bounces < 2u)
                envSample.rgb = uEnvironmentColor;
            vec3 radiance = photon.transmittance * envSample.rgb * uEnvironmentContribution;
            photon.samples++;
            photon.radiance += (radiance - photon.radiance) / float(photon.samples);
            resetPhoton(r, photon);
        } else if (photon.bounces >= uMaxBounces) {
            // max bounces achieved -> only estimate transmittance
            float weightAS = (muAbsorption + muScattering) / uMajorant;
            photon.transmittance *= 1.0 - weightAS;
        } else if (r.y < PAbsorption) {
            // absorption
            float weightA = muAbsorption / (uMajorant * PAbsorption);
            photon.transmittance *= 1.0 - weightA;
        } else if (r.y < PAbsorption + PScattering) {
            // scattering
            r = rand(r);
            float weightS = muScattering / (uMajorant * PScattering);
            photon.transmittance *= volumeSample.rgb * weightS;
            photon.direction = sampleHenyeyGreenstein(uScatteringBias, r, photon.direction);
            photon.bounces++;
        } else {
            // null collision
            float weightN = muNull / (uMajorant * PNull);
            photon.transmittance *= weightN;
            // naveden odboj -> propagiraj foton in določi, ali popolni odboj ali lom
            photon.direction = determineNewPhotonDirection(photon);
        }
        photon.previousDensity = photon.density;
    }

    oPosition = vec4(photon.position, 0);
    oDirection = vec4(photon.direction, float(photon.bounces));
    oTransmittance = vec4(photon.transmittance, 0);
    oRadiance = vec4(photon.radiance, float(photon.samples));
}

// #section MCMRender/vertex

#version 300 es

layout (location = 0) in vec2 aPosition;
out vec2 vPosition;

void main() {
    vPosition = (aPosition + 1.0) * 0.5;
    gl_Position = vec4(aPosition, 0.0, 1.0);
}

// #section MCMRender/fragment

#version 300 es
precision mediump float;

uniform mediump sampler2D uColor;

in vec2 vPosition;
out vec4 oColor;

void main() {
    oColor = vec4(texture(uColor, vPosition).rgb, 1);
}

// #section MCMReset/vertex

#version 300 es

layout (location = 0) in vec2 aPosition;

out vec2 vPosition;

void main() {
    vPosition = aPosition;
    gl_Position = vec4(aPosition, 0.0, 1.0);
}

// #section MCMReset/fragment

#version 300 es
precision mediump float;

#define SQRT3 1.73205080757

@Photon

uniform mat4 uMvpInverseMatrix;
uniform vec2 uInverseResolution;
uniform float uRandSeed;
uniform float uBlur;

uniform vec3 uMinCutPlaneValues;
uniform vec3 uMaxCutPlaneValues;
uniform float uViewCutDistance;

in vec2 vPosition;

layout (location = 0) out vec4 oPosition;
layout (location = 1) out vec4 oDirection;
layout (location = 2) out vec4 oTransmittance;
layout (location = 3) out vec4 oRadiance;

@rand
@unprojectRand
@intersectCube

void main() {
    Photon photon;
    vec3 from, to;
    vec2 randState = rand(vPosition * uRandSeed);
    unprojectRand(randState, vPosition, uMvpInverseMatrix, uInverseResolution, uBlur, from, to);
    photon.direction = normalize(to - from);
    vec2 tbounds = max(limitedIntersectCube(from, photon.direction, uMinCutPlaneValues, uMaxCutPlaneValues), 0.0);
    photon.position = from + (tbounds.x + uViewCutDistance * SQRT3 ) * photon.direction;
    // photon.position = from + tbounds.x * photon.direction;
    photon.transmittance = vec3(1);
    photon.radiance = vec3(1);
    photon.bounces = 0u;
    photon.samples = 0u;
    oPosition = vec4(photon.position, 0);
    oDirection = vec4(photon.direction, float(photon.bounces));
    oTransmittance = vec4(photon.transmittance, 0);
    oRadiance = vec4(photon.radiance, float(photon.samples));
}
