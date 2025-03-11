#ifndef DDGI_PROBE_INDEXING_INCLUDED
#define DDGI_PROBE_INDEXING_INCLUDED

#include "Common/Packing.hlsl"

//------------------------------------------------------------------------
// Probe Indexing Helpers
// Independent of Volume, only affected by coordinate axis type
//------------------------------------------------------------------------

// Get the number of Probes in each Slice along the vertical axis
int DDGIGetProbesPerPlane(int3 probeCounts)
{
#if 1
    // Left or Right Y-UP
    return (probeCounts.x * probeCounts.z);
#elif 0
    // Left or Right Z-UP
    return (probeCounts.x * probeCounts.y);
#endif
}

int DDGIGetPlaneIndex(int3 probeCoords)
{
#if 0
    // Left or Right Z-UP
    return probeCoords.z;
#else
    // Left or Right Y-UP
    return probeCoords.y;
#endif
}

int DDGIGetProbeIndexInPlane(int3 probeCoords, int3 probeCounts)
{
#if 1
    // Left or Right Y-UP
    return probeCoords.x + (probeCounts.x * probeCoords.z);
#elif 0
    // Left Z-UP
    return probeCoords.y + (probeCounts.y * probeCoords.x);
#elif 0
    // Right Z-UP
    return probeCoords.x + (probeCounts.x * probeCoords.y);
#endif
}

//------------------------------------------------------------------------
// Probe Index
// Affected by the current Volume, using Helper functions to get the correct 1D index
//------------------------------------------------------------------------

int DDGIGetProbeIndex(int3 probeCoords)
{
    int probesPerPlane = DDGIGetProbesPerPlane(_ProbeCount);
    int planeIndex = DDGIGetPlaneIndex(probeCoords);
    int probeIndexInPlane = DDGIGetProbeIndexInPlane(probeCoords, _ProbeCount);

    return (planeIndex * probesPerPlane) + probeIndexInPlane;
}

//------------------------------------------------------------------------
// Probe Grid Coordinates
// Affected by the current Volume, using Helper functions to get the correct 3D index
//------------------------------------------------------------------------

int3 DDGIGetProbeCoords(int probeIndex)
{
    int3 probeCoords;

#if 1
    // Left or Right Y-UP
    probeCoords.x = probeIndex % _ProbeCount.x;
    probeCoords.y = probeIndex / (_ProbeCount.x * _ProbeCount.z);
    probeCoords.z = (probeIndex / _ProbeCount.x) % _ProbeCount.z;
#elif 0
    // Left Z-UP
    probeCoords.x = (probeIndex / _ProbeCount.y) % _ProbeCount.x;
    probeCoords.y = probeIndex % _ProbeCount.y;
    probeCoords.z = probeIndex / (_ProbeCount.x * _ProbeCount.y);
#elif 0
    // Right Z-UP
    probeCoords.x = probeIndex % _ProbeCount.x;
    probeCoords.y = (probeIndex / _ProbeCount.x) % _ProbeCount.y;
    probeCoords.z = probeIndex / (_ProbeCount.y * _ProbeCount.x);
#endif

    return probeCoords;
}

//------------------------------------------------------------------------
// Texture Coordinates
//------------------------------------------------------------------------

uint3 DDGIGetRayDataTexelCoords(int rayIndex, int probeIndex)
{
    int probesPerPlane = DDGIGetProbesPerPlane(_ProbeCount);

    uint3 coords;
    coords.x = rayIndex;
    coords.z = probeIndex / probesPerPlane;
    coords.y = probeIndex - (coords.z * probesPerPlane);

    return coords;
}

uint3 DDGIGetProbeTexelCoordsOneByOne(int probeIndex)
{
    // Find the probe's plane index
    int probesPerPlane  = DDGIGetProbesPerPlane(_ProbeCount);
    int planeIndex      = int(probeIndex / probesPerPlane);

#if 1
    // Left or Right Y-UP
    int x = (probeIndex % _ProbeCount.x);
    int y = (probeIndex / _ProbeCount.x) % _ProbeCount.z;
#elif 0
    // Left Z-UP
    int x = (probeIndex % _ProbeCount.y);
    int y = (probeIndex / _ProbeCount.y) % _ProbeCount.x;
#elif 0
    // Right Z-UP
    int x = (probeIndex % _ProbeCount.x);
    int y = (probeIndex / _ProbeCount.x) % _ProbeCount.y;
#endif

    return uint3(x, y, planeIndex);
}
// Used to get the top-left Texel Coord of the current Probe, considering expansion. 
// The returned result is the Border Texel, so you need to add uint3(1,1,0) to correctly index the actual content coordinates.
uint3 DDGIGetProbeBaseTexelCoords(int probeIndex, int numProbeInteriorTexels)
{
    uint3 coords = DDGIGetProbeTexelCoordsOneByOne(probeIndex);

    int numProbeTexels = numProbeInteriorTexels + 2;

    return uint3(coords.x * numProbeTexels, coords.y * numProbeTexels, coords.z);
}

float3 DDGIGetProbeUV(int probeIndex, float2 octantCoordinates, int numProbeInteriorTexels)
{
    // Get the probe's texel coordinates, assuming one texel per probe
    uint3 coords = DDGIGetProbeTexelCoordsOneByOne(probeIndex);

    // Add the border texels to get the total texels per probe
    float numProbeTexels = (numProbeInteriorTexels + 2.f);

#if 1
    // Left or Right Y-UP
    float textureWidth = numProbeTexels * _ProbeCount.x;
    float textureHeight = numProbeTexels * _ProbeCount.z;
#elif 0
    // Left Z-UP
    float textureWidth = numProbeTexels * _ProbeCount.y;
    float textureHeight = numProbeTexels * _ProbeCount.x;
#elif 0
    // Right Z-UP
    float textureWidth = numProbeTexels * _ProbeCount.x;
    float textureHeight = numProbeTexels * _ProbeCount.y;
#endif

    // Move to the center of the probe and move to the octant texel before normalizing
    float2 uv   = float2(coords.x * numProbeTexels, coords.y * numProbeTexels) + (numProbeTexels * 0.5f);
    uv          += octantCoordinates.xy * ((float)numProbeInteriorTexels * 0.5f);
    uv          /= float2(textureWidth, textureHeight);
    
    return float3(uv, coords.z);
}

float3 DDGIGetProbeUV(int probeIndex, float3 direction, int numProbeInteriorTexels)
{
    float2 octantCoordinates = EncodeNormalOctahedron(normalize(direction));
    return DDGIGetProbeUV(probeIndex, octantCoordinates, numProbeInteriorTexels);
}

//------------------------------------------------------------------------
// Probe Classification
//------------------------------------------------------------------------



//------------------------------------------------------------------------
// Infinite Scrolling
//------------------------------------------------------------------------

int DDGIGetScrollingProbeIndex(int3 probeCoords)
{
    // TODO: Support Scroll
    const float probeScrollOffsets = 0.0f;
    
    return DDGIGetProbeIndex(((probeCoords + probeScrollOffsets + _ProbeCount) % _ProbeCount));
}

#endif