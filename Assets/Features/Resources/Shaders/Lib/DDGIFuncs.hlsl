#ifndef DDGI_FUNCS
#define DDGI_FUNCS

#include "Common/Packing.hlsl"

//------------------------------------------------------------------------
// Ray tracing payload helper.
//------------------------------------------------------------------------
DDGIPayload GetPrimaryPayload()
{
	const DDGIPayload payload = (DDGIPayload) 0;

	return payload;
}

DDGIPayload GetShadowPayload()
{
	DDGIPayload payload = (DDGIPayload) 0;

	payload.isShadowPayload = true;
	payload.isInShadow		= false;

	return payload;
}


//------------------------------------------------------------------------
// Math Utility
//------------------------------------------------------------------------

float Min(float2 v) { return min(v.x, v.y); }
float Min(float3 v) { return min(Min(v.xy), v.z); }
float Min(float4 v) { return min(Min(v.xyz), v.w);}
float Max(float2 v) { return max(v.x, v.y); }
float Max(float3 v) { return max(Max(v.xy), v.z); }
float Max(float4 v) { return max(Max(v.xyz), v.w);}

float Pow2(float x) { return x * x; }
float Pow3(float x) { return x * x * x; }

// Outputs the rotated direction (in radians)
// Reference: Unity Shader Graph
float3 RotateAboutAxisInRadians(float3 In, float3 Axis, float Rotation)
{
    float s = sin(Rotation);
    float c = cos(Rotation);
    float one_minus_c = 1.0 - c;

    Axis = normalize(Axis);
    float3x3 rot_mat =
    { one_minus_c * Axis.x * Axis.x + c, one_minus_c * Axis.x * Axis.y - Axis.z * s, one_minus_c * Axis.z * Axis.x + Axis.y * s,
        one_minus_c * Axis.x * Axis.y + Axis.z * s, one_minus_c * Axis.y * Axis.y + c, one_minus_c * Axis.y * Axis.z - Axis.x * s,
        one_minus_c * Axis.z * Axis.x - Axis.y * s, one_minus_c * Axis.y * Axis.z + Axis.x * s, one_minus_c * Axis.z * Axis.z + c
    };
    return mul(rot_mat, In);
}

// Outputs the rotation matrix (in radians)
float3x3 AngleAxis3x3(float angle, float3 axis)
{
	// Rotation with angle (in radians) and axis
	float c, s;
	sincos(angle, s, c);

	float t = 1 - c;
	float x = axis.x;
	float y = axis.y;
	float z = axis.z;

	return float3x3(
		t * x * x + c, t * x * y - s * z, t * x * z + s * y,
		t * x * y + s * z, t * y * y + c, t * y * z - s * x,
		t * x * z - s * y, t * y * z + s * x, t * z * z + c
		);
}


//------------------------------------------------------------------------
// Quaternion Helpers
//------------------------------------------------------------------------

/**
 * Rotate vector v with quaternion q.
 */
float3 DDGIQuaternionRotate(float3 v, float4 q)
{
	float3 b = q.xyz;
	float b2 = dot(b, b);
	return (v * (q.w * q.w - b2) + b * (dot(v, b) * 2.f) + cross(b, v) * (q.w * 2.f));
}

/**
 * Quaternion conjugate.
 * For unit quaternions, conjugate equals inverse.
 * Use this to create a quaternion that rotates in the opposite direction.
 */
float4 DDGIQuaternionConjugate(float4 q)
{
	return float4(-q.xyz, q.w);
}


//------------------------------------------------------------------------
// Randomize Functions
//------------------------------------------------------------------------

// Ray Tracing Gems 2: Essential Ray Generation Shaders
float3 SphericalFibonacci(float i, float n)
{
	const float PHI = sqrt(5) * 0.5f + 0.5f;
	float fraction	= (i * (PHI - 1)) - floor(i * (PHI - 1));
	float phi		= 2.0f * PI * fraction;
	float cosTheta	= 1.0f - (2.0f * i + 1.0f) * (1.0f / n);
	float sinTheta	= sqrt(saturate(1.0 - cosTheta * cosTheta));
	
	return float3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
}

float3 DDGIGetProbeRayDirection(int rayIndex)
{
	bool isFixedRay = false;
	int sampleIndex = rayIndex;
	int numRays		= _RaysPerProbe;
	
	if ((DDGI_PROBE_RELOCATION == DDGI_PROBE_RELOCATION_ON) || (DDGI_PROBE_CLASSIFICATION == DDGI_PROBE_CLASSIFICATION_ON))
	{
		isFixedRay  = (rayIndex < RTXGI_DDGI_NUM_FIXED_RAYS);
		sampleIndex = isFixedRay ? rayIndex : (rayIndex - RTXGI_DDGI_NUM_FIXED_RAYS);
		numRays		= isFixedRay ? RTXGI_DDGI_NUM_FIXED_RAYS : (numRays - RTXGI_DDGI_NUM_FIXED_RAYS);
	}

	// Get a ray direction on the sphere
	float3 direction = SphericalFibonacci(sampleIndex, numRays);

	// Don't rotate fixed rays so relocation/classification are temporally stable
	if (isFixedRay) return normalize(direction);

	// Apply Rotation
	float3 randomDirection = RotateAboutAxisInRadians(direction, _RandomVector, _RandomAngle);
	return normalize(randomDirection);
}


//------------------------------------------------------------------------
// Light Fetcher
//------------------------------------------------------------------------

Light GetDDGIDirectionalLight(int index)
{
	DirectionalLight directionalLight = DirectionalLightBuffer[index];

	Light light;
	light.direction				= directionalLight.direction.xyz;
	light.color					= directionalLight.color.rgb;
	light.distanceAttenuation	= 1.0f;
	light.shadowAttenuation		= 1.0f;
	light.layerMask				= 0;

	return light;
}

Light GetDDGIPunctualLight(int index, float3 positionWS)
{
	PunctualLight punctualLight = PunctualLightBuffer[index];
	float4 lightPositionWS = punctualLight.position;
	float3 color = punctualLight.color.rgb;
	float4 distanceAndSpotAttenuation = punctualLight.distanceAndSpotAttenuation;
	float4 spotDirection = punctualLight.spotDirection;
	
	float3 lightVector	= lightPositionWS.xyz - positionWS * lightPositionWS.w;
	float distanceSqr	= max(dot(lightVector, lightVector), HALF_MIN);

	half3 lightDirection = half3(lightVector * rsqrt(distanceSqr));
	float attenuation	 = DistanceAttenuation(distanceSqr, distanceAndSpotAttenuation.xy) * AngleAttenuation(spotDirection.xyz, lightDirection, distanceAndSpotAttenuation.zw);

    // We use ray tracing to determine shadows, here shadowAttenuation is set to 1.
	Light light;
	light.direction				= lightDirection;
	light.distanceAttenuation	= attenuation;
	light.shadowAttenuation		= 1.0;
	light.color					= color.rgb;
	light.layerMask				= 0;

	return light;
}


//------------------------------------------------------------------------
// Probe Data Fetcher
//------------------------------------------------------------------------

#if defined(DDGI_VISUALIZATION) || defined(DDGI_RAYTRACING) || defined(FORWARD_USE_DDGI)
	float3 DDGILoadProbeDataOffset(uint3 coords)
	{
		return LOAD_TEXTURE2D_ARRAY_LOD(_ProbeData, coords.xy, coords.z, 0).xyz * _ProbeSize;	
	}

	int DDGILoadProbeState(uint3 coords)
	{
	    int state = DDGI_PROBE_STATE_ACTIVE;
		if(DDGI_PROBE_CLASSIFICATION == DDGI_PROBE_CLASSIFICATION_ON)
		{
			state = (int)LOAD_TEXTURE2D_ARRAY_LOD(_ProbeData, coords.xy, coords.z, 0).a;
		} 

		return state;
	}
#else
	// We use Texture2DArray in visualization, Texture2DArray dont support these function.
	float3 DDGILoadProbeDataOffset(uint3 coords)
	{
		return _ProbeData[coords].xyz * _ProbeSize;
	}

	void DDGIStoreProbeDataOffset(uint3 coords, float3 wsOffset)
	{
		// A-Component is useless now.
		_ProbeData[coords] = float4(wsOffset / _ProbeSize, 1.0f);
	}

	int DDGILoadProbeState(uint3 coords)
	{
		int state = DDGI_PROBE_STATE_ACTIVE;
		if(DDGI_PROBE_CLASSIFICATION == DDGI_PROBE_CLASSIFICATION_ON)
		{
			state = (int)_ProbeData[coords].a;
		}

		return state;
	}
#endif


//------------------------------------------------------------------------
// Probe World Position
//------------------------------------------------------------------------

// Get the world position of a probe based on its 3D grid coordinates
float3 DDGIGetProbeWorldPosition(uint3 gridCoord)
{
    // float3 probeWorldPosition = _StartPosition + _ProbeSize * gridCoord; // No Rotation Implementation

    const float3 probeSpaceWorldPosition = gridCoord * _ProbeSize;
    const float3 probeVolumeExtents = (_ProbeSize * (_ProbeCount - 1)) * 0.5f; // Actually ddgiVolumeCpu - extents

    // Rotate Probe
    float3 probeWorldPosition = probeSpaceWorldPosition - probeVolumeExtents; // Convert [0,n] coordinates to [-n/2, n/2] for center rotation
    probeWorldPosition = DDGIQuaternionRotate(probeWorldPosition, _ProbeRotation) + probeVolumeExtents;
    probeWorldPosition += _StartPosition;

    // Ray tracing shaders will use this function. According to the link below, ray tracing shader branches are still in the plan, 
    // meaning we cannot use variants, so use a variable to check if it's enabled or not.
    // https://portal.productboard.com/unity/1-unity-platform-rendering-visual-effects/tabs/125-shader-system
    if (DDGI_PROBE_RELOCATION == DDGI_PROBE_RELOCATION_ON)
    {
        // Since when sampling tex2DArray, the z-component of the sample coordinate actually corresponds to the y-component of gridCoord, 
        // an extra inversion is needed here.
        int probeIndex = DDGIGetProbeIndex(gridCoord);
        uint3 probeDataTexelCoord = DDGIGetProbeTexelCoordsOneByOne(probeIndex);
        probeWorldPosition += DDGILoadProbeDataOffset(probeDataTexelCoord);
    }

    return probeWorldPosition;
}

// Get the world position of a probe based on its 1D grid index
float3 DDGIGetProbeWorldPosition(uint probeIndex)
{
    uint3 gridCoord = DDGIGetProbeCoords(probeIndex);
    return DDGIGetProbeWorldPosition(gridCoord);
}

// Given a world space position P, return the base grid coordinates of the probe related to this position 
// (used to determine the probe grid block that P belongs to)
uint3 DDGIGetBaseGridCoords(float3 worldPos)
{
	const float3 probeVolumeExtents  = (_ProbeSize * (_ProbeCount - 1)) * 0.5f;
	const float3 probeVolumeCenter   = _StartPosition + probeVolumeExtents;

	float3 position = worldPos - probeVolumeCenter;
	position = DDGIQuaternionRotate(position, DDGIQuaternionConjugate(_ProbeRotation));
	position += probeVolumeExtents; // Transform to [0,n]

	uint3 probeCoords = uint3(position / _ProbeSize);
	probeCoords = clamp(probeCoords, uint3(0,0,0), uint3(_ProbeCount) - uint3(1, 1, 1));

	return probeCoords;
	
	// return clamp(uint3((worldPos - _StartPosition) / _ProbeSize), uint3(0, 0, 0), uint3(_ProbeCount) - uint3(1, 1, 1)); // No Rotation Implementation.
}

//------------------------------------------------------------------------
// Runtime Probe Sampling
//------------------------------------------------------------------------

float3 ComputeBias(float3 normal, float3 viewDirection, float b = 0.2f)
{
	#if 0
		// Arida Implementation.
	    const float normalBiasMultiplier = 0.2f;
	    const float viewBiasMultiplier = 0.8f;
	    const float axialDistanceMultiplier = 0.75f;
	    return (normal * normalBiasMultiplier + viewDirection * viewBiasMultiplier) * axialDistanceMultiplier * Min(_ProbeSize) * b;
	#else
		// NVIDIA Implementation.
		return (normal * _NormalBiasMultiplier + viewDirection * _ViewBiasMultiplier);
	#endif
}

float DDGIGetVolumeBlendWeight(float3 worldPosition)
{
	const float3 probeVolumeExtents  = (_ProbeSize * (_ProbeCount - 1)) * 0.5f;
	const float3 probeVolumeCenter   = _StartPosition + probeVolumeExtents;

	float3 position = worldPosition - probeVolumeCenter;
	position = abs(DDGIQuaternionRotate(position, DDGIQuaternionConjugate(_ProbeRotation)));

	float3 delta = position - probeVolumeExtents;
	if(all(delta < 0.0f)) return 1.0f;

	float volumeBlendWeight = 1.0f;
	volumeBlendWeight *= (1.0f - saturate(delta.x / _ProbeSize.x));
	volumeBlendWeight *= (1.0f - saturate(delta.y / _ProbeSize.y));
	volumeBlendWeight *= (1.0f - saturate(delta.z / _ProbeSize.z));

	return volumeBlendWeight;
}

//https://github.com/simco50/D3D12_Research/blob/master/D3D12/Resources/Shaders/RayTracing/DDGICommon.hlsli
float3 SampleDDGIIrradiance(float3 P, float3 N, float3 Wo)
{
	float3 direction		= N;
	float3 biasedPosition	= P;
	float3 unbiasedPosition = P;
	float  volumeWeight		= 1.0f;

	biasedPosition += ComputeBias(direction, -Wo);

    // When the shading point is outside the Volume area, we return early.
    // When the shading point approaches the Volume boundary (but doesn't exceed the volume area), 
    // we perform a smooth transition for its irradiance.
    volumeWeight = DDGIGetVolumeBlendWeight(biasedPosition);
    if (volumeWeight <= 0.0f) return 0.0f;

    // When computing relativeCoordinates, the position needs to be offset (reference to NVIDIA).
    // If the position is offset here (Adria's implementation), it will cause mesh artifacts during trilinear interpolation.
    // position += ComputeBias(direction, -Wo);

	const uint3  baseProbeCoords	= DDGIGetBaseGridCoords(biasedPosition);
	const float3 baseProbePosition	= DDGIGetProbeWorldPosition(baseProbeCoords);

	float3 gridSpaceDistance = biasedPosition - baseProbePosition;
	gridSpaceDistance		 = DDGIQuaternionRotate(gridSpaceDistance, DDGIQuaternionConjugate(_ProbeRotation));
	const float3 alpha		 = saturate(gridSpaceDistance / _ProbeSize);

	float3 sumIrradiance = 0;
	float  sumWeight	 = 0;

	for (uint j = 0; j < 8; ++j)
	{
		const uint3 indexOffset = uint3(j, j >> 1u, j >> 2u) & 1u;

		const uint3 probeCoords		= clamp(baseProbeCoords + indexOffset, 0, _ProbeCount - 1);
		const float3 probePosition	= DDGIGetProbeWorldPosition(probeCoords);
		const uint probeIndex		= DDGIGetProbeIndex(probeCoords);

        // If the probe is not enabled, do not interpolate.
        const uint3 probeDataCoords = DDGIGetProbeTexelCoordsOneByOne(probeIndex);
		const int   probeState      = DDGILoadProbeState(probeDataCoords);
		if(probeState == DDGI_PROBE_STATE_INACTIVE) continue;

		float3 relativeProbePosition = biasedPosition - probePosition;
		float3 probeDirection		 = -normalize(relativeProbePosition);

		float3 trilinear		= max(0.001f, lerp(1.0f - alpha, alpha, indexOffset));
        float trilinearWeight	= (trilinear.x * trilinear.y * trilinear.z);

		float weight = 1.0f;

		// --------------------------------
		// Backface Weighting
		// --------------------------------
		#if 0
			// Adria Implementation.
			weight *= saturate(dot(probeDirection, direction));
		#else
			// NVIDIA Implementation.
			const float wrapShading = dot(normalize(probePosition - unbiasedPosition), direction) * 0.5f + 0.5f;
			weight *= (wrapShading * wrapShading) + 0.2f;
		#endif

		// --------------------------------
		// Chebyshev Visibility Test
		// --------------------------------
		float3 probeDistanceUV	= DDGIGetProbeUV(probeIndex, -probeDirection, PROBE_DISTANCE_TEXELS);
		float  probeDistance	= length(relativeProbePosition);
		// https://developer.download.nvidia.com/SDK/10/direct3d/Source/VarianceShadowMapping/Doc/VarianceShadowMapping.pdf
		float2 moments = SAMPLE_TEXTURE2D_ARRAY_LOD(_ProbeDistanceHistory, sampler_LinearClamp, probeDistanceUV.xy, probeDistanceUV.z, 0).xy;
		float variance = abs(Pow2(moments.x) - moments.y);
		float chebyshev = 1.0f;
		if(probeDistance > moments.x)
		{
			float mD = moments.x - probeDistance;
			chebyshev = variance / (variance + Pow2(mD));
			chebyshev = max(Pow3(chebyshev), 0.0);
		}
		weight *= max(chebyshev, 0.05f);
		
		weight = max(0.000001f, weight);

		// --------------------------------
		// Threshold and Trilinear Weight.
		// --------------------------------
		const float crushThreshold = 0.2f;
		if (weight < crushThreshold)
		{
			weight *= weight * weight * (1.0f / Pow2(crushThreshold));
		}
		weight *= trilinearWeight;

		float3 probeIrradianceUV = DDGIGetProbeUV(probeIndex, direction, PROBE_IRRADIANCE_TEXELS);
		float3 irradiance		 = SAMPLE_TEXTURE2D_ARRAY_LOD(_ProbeIrradianceHistory, sampler_LinearClamp, probeIrradianceUV.xy, probeIrradianceUV.z, 0).rgb;
		irradiance				 = pow(irradiance, 2.5f); // Gamma Correct.

		sumIrradiance += irradiance * weight;
		sumWeight	  += weight;
	}
	
	if(sumWeight == 0) return 0.0f;

	sumIrradiance *= (1.0f / sumWeight);
	sumIrradiance *= sumIrradiance;
	sumIrradiance *= DDGI_2PI;
	sumIrradiance *= _IndirectIntensity;
	
	return sumIrradiance * volumeWeight;
}


// ---------------------------------------
// Legacy (Ignored)
// ---------------------------------------

/*// Get the world position of a probe based on its 3D grid coordinates (considering Relocation)
float3 DDGIGetRelocatedProbeWorldPosition(int3 probeCoords)
{
    float3 probeWorldPosition = DDGIGetProbeWorldPosition(probeCoords);

    int probeIndex        = DDGIGetProbeIndex(probeCoords);
    uint3 coords          = DDGIGetProbeTexelCoordsOneByOne(probeIndex);
    probeWorldPosition   += DDGILoadProbeDataOffset(coords);

    return probeWorldPosition;
}

// Get the world position of a probe based on its 1D grid index (considering Relocation)
float3 DDGIGetRelocatedProbeWorldPosition(int probeIndex)
{
    int3 probeCoords = DDGIGetProbeCoords(probeIndex);
    return DDGIGetRelocatedProbeWorldPosition(probeCoords);
}*/

#endif