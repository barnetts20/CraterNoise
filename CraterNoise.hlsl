struct CraterNoise {
    struct CraterCurve{
        float X;
        float K1;
        float K2;        
        float floorMulti;
        float floorExp;
        float floorOffset;
        float cavityMulti;
        float cavityExp;
        float ridgeMulti;
        float ridgeExp;
        float finalMulti;
        
        float Smooth(float a, float b, float k) {
            return -log(exp(-k * a) + exp(-k * b)) / k;
        }
        
        float CraterFloor() {
            float f1 = Smooth(1 - pow(abs(4 * X - 4), floorExp), 0, -K2);
            return floorMulti * Smooth(f1, 0.0, K2) + floorOffset;
        }

        float CraterCavity() {
            return (X >= 0 && X <= 2) ? -cavityMulti * pow(0.5 * (sin((X + 1.5) * PI) + 1), cavityExp) + 2.0 : 2.0;
        }

        float CraterRidge() {
            return (X >= 0 && X <= 2) ? ridgeMulti * pow(0.5 * (sin((X + 1.5) * PI) + 1), ridgeExp) : 0.0;
        }

        float Crater() {
            float cavityFloor = Smooth(CraterFloor(), CraterCavity(), -K1);
            return (Smooth(CraterRidge(), cavityFloor, K1) + max(0.0, (2 - K1))) * finalMulti;
        }

        void SeedCrater(float3 randomFactor){
            // Example of randomizing the parameters based on randomFactor input
            randomFactor = frac(randomFactor);
            floorMulti = lerp(-4, -1, randomFactor.x);
            floorExp = lerp(-.2, .05, randomFactor.y);
            floorOffset = lerp(-2, -0.1, randomFactor.z);
            cavityMulti = lerp(2, 4, randomFactor.x);
            cavityExp = lerp(1.0, 4.0, randomFactor.y);
            ridgeMulti = lerp(2.0, 4.0, randomFactor.z);
            ridgeExp = lerp(1.0, 4.0, randomFactor.x);
            finalMulti = lerp(0.25, 1.75, randomFactor.y);
        }

        void SeedVolcano(float3 randomFactor){
            randomFactor = frac(randomFactor);
            floorMulti = lerp(-4, -1, randomFactor.x);
            floorExp = lerp(-.2, -.05, randomFactor.y);
            floorOffset = lerp(-10, -8.0, randomFactor.z);
            cavityMulti = lerp(4, 8, randomFactor.x);
            cavityExp = lerp(5.0, 10.0, randomFactor.y);
            ridgeMulti = lerp(3.0, 6.0, randomFactor.z);
            ridgeExp = lerp(1.0, 3.0, randomFactor.x);
            finalMulti = lerp(0.5, 1.5, randomFactor.y);
        }
    };

    int OCTAVES;
    float OCTAVE_FREQUENCY_SCALE;
    float OCTAVE_AMPLITUDE_SCALE;

    float3 Noise(float3 p) {
        return frac(sin(float3(dot(p, float3(127.1, 311.7, 591.1)), 
                                dot(p, float3(269.5, 183.3, 113.5)), 
                                dot(p, float3(419.2, 371.9, 4297.7)))) * 43758.5454453);
    }
    
    // This function calculates the intersecting circle between two overlapping spheres
    float4 CalculateTwoSphereIntersect(float3 sphereCenter1, float sphereRadius1, float3 sphereCenter2, float sphereRadius2) {
        // Calculate the distance between the centers of the two spheres
        float distanceBetweenCenters = length(sphereCenter2 - sphereCenter1);
        // Check if the spheres do not intersect or are completely contained within each other
        if (distanceBetweenCenters > sphereRadius1 + sphereRadius2 || distanceBetweenCenters < abs(sphereRadius1 - sphereRadius2)) { return float4(sphereCenter2, 0); } // No intersection
        // Calculate the radius of the intersection circle
        float intersectionRadius = sqrt(sphereRadius1 * sphereRadius1 - pow((distanceBetweenCenters * distanceBetweenCenters + sphereRadius1 * sphereRadius1 - sphereRadius2 * sphereRadius2) / (2.0 * distanceBetweenCenters), 2.0));
        // Calculate the center of the intersection circle
        float3 intersectionCenter = sphereCenter1 + (sphereCenter2 - sphereCenter1) * (sphereRadius1 * sphereRadius1 - sphereRadius2 * sphereRadius2 + distanceBetweenCenters * distanceBetweenCenters) / (2.0 * distanceBetweenCenters * distanceBetweenCenters);
        return float4(intersectionCenter, max(0.0, intersectionRadius));
    }
    
    // in pos should be a scaled value that starts with the unit sphere
    // This function generates a 3D voronoi grid, 
    // and for each cell, creates a sphere starting at its center with a radius of the distance of the center to closest cell wall. 
    // This garauntees no overlaps.
    // Next it generates a 2 sphere circle overlap of the voronoi sphere and the sphere formed by center 0,0,0 and radius equal to the length of the sample position
    // It then uses these values to generate displacement, ejecta, distance, and max distance values in rbgw 
    // These layers are good for creating features like craters or volcanos.
    float4 VoronoiCraterNoise(in float3 samplePos) {
        float3 integerPos = floor(samplePos);
        float3 fracPos = frac(samplePos);
        float3 closestFeatureCenter = float3(0.0, 0.0, 0.0);
        float3 closestFeaturePoint = float3(0.0,0.0,0.0);
        float3 closestDiff = float3(0.0,0.0,0.0);;
        float minEdgeDistance = 8.0;
        float minGradientDistance = 8.0;
        for (int k = -1; k <= 1; k++) {
            for (int j = -1; j <= 1; j++) {
                for (int i = -1; i <= 1; i++) {
                    float3 neighborCell = float3(i, j, k);
                    float3 featurePoint = neighborCell + Noise(integerPos + neighborCell);
                    float3 diff = featurePoint - fracPos;
                    float distance = dot(diff, diff);
                    
                    if (distance < minEdgeDistance) {
                        minEdgeDistance = distance;
                        minGradientDistance = distance;
                        closestDiff = diff;
                        closestFeatureCenter = featurePoint;
                        closestFeaturePoint = neighborCell;
                    }
                }
            }
        }

        minEdgeDistance = 8.0;
        minGradientDistance = 8.0;
        for (int k = -2; k <= 2; k++) {
            for (int j = -2; j <= 2; j++) {
                for (int i = -2; i <= 2; i++) {
                    float3 neighborCell = float3(i, j, k) + closestFeaturePoint;
                    float3 neighborFeaturePoint = neighborCell + Noise(integerPos + neighborCell);
                    float distToNeighbor = distance(neighborFeaturePoint, closestFeatureCenter);
                    float3 diff = neighborFeaturePoint - fracPos;
                    if (distToNeighbor > 0.0) {
                        minEdgeDistance = min(minEdgeDistance, distToNeighbor);
                        minGradientDistance = min(minGradientDistance, dot(0.5 * (closestDiff + diff), normalize(diff - closestDiff)));
                    }
                }
            }
        }
        float voronoiSphereRadius = minEdgeDistance * 0.5;
        float3 voronoiCellCenter = closestFeatureCenter + integerPos;
        float sampleSphereRadius = length(samplePos);
        float3 sampleSphereCenter = float3(0, 0, 0);
        // Calculate the intersection sphere using the CalculateTwoSphereIntersect function
        float4 intersectionSphere = CalculateTwoSphereIntersect(sampleSphereCenter, sampleSphereRadius, voronoiCellCenter, voronoiSphereRadius);
        // Calculate displacement, ejecta, sphere distance, and max distance
        float vSphereDistance = max(0, 1 - (distance(intersectionSphere.xyz, samplePos) / intersectionSphere.w));
        CraterCurve cc;
        cc.X = saturate(vSphereDistance);
        cc.K1 = 5;
        cc.K2 = 5;
        cc.SeedCrater(closestFeaturePoint);
        float displacement = cc.Crater();

        float ejecta = CalculateEjecta(intersectionSphere.xyz, intersectionSphere.w, float3(0, 0, 0), samplePos);
        float maxDistance = intersectionSphere.w;
        return float4(displacement, ejecta, vSphereDistance, minGradientDistance);
    }

    //Procedural Ejecta function
    float CalculateEjecta(float3 gradientCenter, float gradientRadius, float3 objectCenter, float3 samplePosition) {
        // Calculate the surface normal at the gradient center
        float3 N = normalize(gradientCenter - objectCenter);
        // Calculate tangent and bitangent vectors to define the tangent plane
        float3 up = abs(N.y) < 0.999 ? float3(0.0, 1.0, 0.0) : float3(1.0, 0.0, 0.0);
        float3 tangent = normalize(cross(up, N));
        float3 bitangent = cross(N, tangent);
        // Project the sample position onto the tangent plane
        float3 toSample = samplePosition - gradientCenter;
        float u = dot(toSample, tangent);
        float v = dot(toSample, bitangent);
        // Radial distance from the center of the gradient sphere
        float distanceFromCenter = length(float2(u, v));
        // Normalize the distance in the range [0, 1]
        float normalizedDistance = distanceFromCenter / gradientRadius;
        // Azimuthal angle (angle around the normal vector in the tangent plane)
        float angle = degrees(atan2(v, u));
        angle = fmod(angle + 360.0, 360.0); // Ensure angle is in the range [0, 360)
        // Ejecta strength falls off with distance
        float ejectaStrength = exp(-normalizedDistance * 3.0);  // Example: Exponential falloff
        // Generate a random factor based on the integer part of the gradient center to vary properties per ejecta site
        float3 intGradientCenter = floor(gradientCenter);
        float randomFactor = Noise(intGradientCenter).x;  // Use the x component of the noise result as a random float
        randomFactor = frac(randomFactor * 10.0);  // Keep it in the range [0, 1]
        // Determine the number of angular segments
        int numSegments = lerp(5, 10, randomFactor); // Random between 5 and 10 mirrored segments
        // Initialize a variable to check if the angle falls within any segment
        bool angleInSegment = false;
        float segmentStrength = 1.0;
        float segmentFrequency = 1.0;

        for (int i = 0; i < numSegments; i++) {
            // Generate a random starting angle for the segment, ensuring it covers the full 360 degrees
            float segmentStart = frac(sin(dot(intGradientCenter + float3(i, i, i), float3(12.9898, 78.233, 45.164))) * 43758.5453123) * 360.0;
            segmentStart = fmod(segmentStart, 360.0); // Ensure it stays within [0, 360)
            // Generate a random segment width between 15 to 25 degrees
            float segmentWidth = lerp(5.0, 25.0, frac(sin(dot(intGradientCenter + float3(i + 1, i + 1, i + 1), float3(12.9898, 78.233, 45.164))) * 43758.5453123));
            // Calculate the mirrored start angle
            float mirroredStart = fmod(segmentStart + 180.0, 360.0);
            // Check if the current angle falls within this segment or its mirrored segment
            if ((angle >= segmentStart && angle <= segmentStart + segmentWidth) ||
                (angle >= mirroredStart && angle <= mirroredStart + segmentWidth)) {
                angleInSegment = true;
                // Increase strength variance within segments
                segmentStrength = lerp(0.25, 4.0, frac(sin(dot(intGradientCenter + float3(i + 2, i + 2, i + 2), float3(12.9898, 78.233, 45.164))) * 43758.5453123));
                // Bind frequency randomization to the segment
                segmentFrequency = lerp(.10, 4.0, frac(sin(dot(intGradientCenter + float3(i + 3, i + 3, i + 3), float3(12.9898, 78.233, 45.164))) * 43758.5453123));
                break;
            }
        }
        // If the angle falls within any of the segments, apply the segment's strength and frequency
        if (angleInSegment) {
            ejectaStrength *= segmentStrength;
            // Calculate the streak pattern with segment-specific variations
            float streakPattern = smoothstep(0.5 - segmentFrequency * 0.1, 0.5 + segmentFrequency * 0.1, 0.5 + 0.5 * sin(angle * segmentFrequency));
            // Combine streaks with the radial falloff and segment-specific strength
            float ejecta = ejectaStrength * streakPattern;
            return saturate(ejecta);
        }
        return 0.0;
    }

    //Finally, fractal brownian motion to allow us to generate multiple layers of craters at different scales and intensities
    float4 CraterFBM(float3 pos){
        float3 curSamplePos = pos;
        float curAmplitudeCoeff = 1.0;
        float4 accumulatedCraterData = float4(0.0,0.0,0.0,0.0);    
        for(int i = 0; i < OCTAVES; i++){
            accumulatedCraterData += VoronoiCraterNoise(curSamplePos) * curAmplitudeCoeff; 
            curAmplitudeCoeff *= OCTAVE_AMPLITUDE_SCALE;
            curSamplePos *= OCTAVE_FREQUENCY_SCALE;
        }        
        return accumulatedCraterData;
    }
};

CraterNoise cn;
cn.OCTAVES = Octaves;
cn.OCTAVE_FREQUENCY_SCALE = OctaveFrequencyScale;
cn.OCTAVE_AMPLITUDE_SCALE = OctaveAmplitudeScale;
float3 scaledPosition = WorldPosition/Radius * Frequency;

return cn.CraterFBM(scaledPosition);