using Revise
using FileIO
using GeometryBasics
using BenchmarkTools
using GLMakie

include("Bresenham.jl")

planeNormal = Point3{Float32}(1.0,0.0,0.0)
function getLines(ds::Number, verts::Array{T,1},tris::Vector{Vector{Int64}}) where T<:AbstractPoint
    pointResult = Set{T}()
    lines = Vector{T}()

    mins = [Inf32,Inf32,Inf32]
    maxs = [-Inf32,-Inf32,-Inf32]

    #Search all triangles
    for t in 1:size(tris,1)
        @views getVertExtrema!(verts[tris[t]],maxs,mins)
        minX = round(Int64, mins[1]/ds)
        maxX = round(Int64, maxs[1]/ds)
        @assert(minX<=maxX, "Error in tri bounding boxes")
        #Go from the lowest point in the tri to highest
        if minX != maxX
            for i in minX:maxX
                empty!(pointResult)
                #Each of the three lines
                segmentPlaneIntersection!(verts[tris[t][1]],verts[tris[t][2]],planeNormal,-i*ds,ds, pointResult)
                segmentPlaneIntersection!(verts[tris[t][2]],verts[tris[t][3]],planeNormal,-i*ds,ds, pointResult)
                segmentPlaneIntersection!(verts[tris[t][3]],verts[tris[t][1]],planeNormal,-i*ds,ds, pointResult)

                if length(pointResult)==2
                    append!(lines,pointResult)
                end
            end
        else
            push!(lines, verts[tris[t][1]])
            push!(lines, verts[tris[t][2]])

            push!(lines, verts[tris[t][2]])
            push!(lines, verts[tris[t][3]])

            push!(lines, verts[tris[t][3]])
            push!(lines, verts[tris[t][1]])
        end
    end
    return lines
end

#https://stackoverflow.com/questions/3142469/determining-the-intersection-of-a-triangle-and-a-plane#
function distFromPlane(point::AbstractPoint, planeN::AbstractPoint, planeD::Number)
    return sum(planeN.*point)+planeD
end
function segmentPlaneIntersection!(p1::T,p2::T, planeN::T, planeD::Number,ds::Number, intersection::Set{T}) where T<: AbstractPoint
    d1 = distFromPlane(p1,planeN,planeD)
    d2 = distFromPlane(p2,planeN,planeD)
    p1OnPlane = abs(d1) < eps(Float32) 
    p2OnPlane = abs(d2) < eps(Float32) 
    
    if p1OnPlane 
        push!(intersection,p1)
    end
    if p2OnPlane
        push!(intersection,p2)
    end
    if p1OnPlane && p2OnPlane #Both points intersect
        return
    end
    if d1*d2 > eps(Float32) #Same side of plane
        return
    end
    t = d1/(d1-d2)
    push!(intersection,p1 + t * (p2 - p1))


end

function getVertExtrema(verts::Array{Point3{Float32},1})
    mins = [Inf32,Inf32,Inf32]
    maxs = [-Inf32,-Inf32,-Inf32]
    for v in verts
        for i in 1:3
            mins[i] = min(mins[i],v[i])
            maxs[i] = max(maxs[i],v[i])
        end
    end
    return (mins,maxs)
end
function getVertExtrema!(verts::SubArray{Point3{Float32},1}, maxs::AbstractVector,mins::AbstractVector)
    mins .= Inf32
    maxs .= -Inf32
    for v in verts
        for i in 1:3
            mins[i] = min(mins[i],v[i])
            maxs[i] = max(maxs[i],v[i])
        end
    end
    return (mins,maxs)
end
function getVertExtrema!(verts::Vector{Point3{Float32}}, maxs::AbstractVector,mins::AbstractVector)
    mins .= Inf32
    maxs .= -Inf32
    for v in verts
        for i in 1:3
            mins[i] = min(mins[i],v[i])
            maxs[i] = max(maxs[i],v[i])
        end
    end
    return (mins,maxs)
end

function scanline(arr::Vector{Bool})
    len = length(arr)
    start = -1
    stop = len
    draw = true
    for i in 1:len
        if arr[i]
            start = i
            break
        end
    end
    #No starting pos found, break
    if start == -1
        return arr
    end
    while start < len
        for i in start+1:len
            if arr[i]
                stop = i
                break
            end
        end
        if stop == start
            break
        end
        if draw
            arr[start:stop] .= true
        end 
        start = stop

        draw = !draw
    end
    return arr
end

function voxelMeshSize(bound::AbstractArray, resolution::Integer)
    longest = findmax(bound)[2]
    ppu = resolution/bound[longest]
    return (round.(Int,bound*ppu),longest,1/ppu)
end
function bresenhamf!(y::Integer,z::Integer,args)
    x=args[1]
    grid = args[2]
    grid[x,y,z]=true
end
function voxelizeMesh(mesh::GeometryBasics.AbstractMesh; resolution::Integer=100)
    verts = Array{Point3{Float32},1}()
    normals = Array{Point3{Float32},1}()

    for t in mesh
        for v in t
            push!(verts,v)
            push!(normals,v.normals)
        end
    end
    #Preallocate vectors 
    mins = [Inf32,Inf32,Inf32]
    maxs = [-Inf32,-Inf32,-Inf32]

    numTris = size(mesh,1)
    tris = [[3*i-2,3*i-1,3*i] for i in 1:numTris]

    #Find the extrema for the mesh
    mins,maxs=getVertExtrema!(verts,maxs,mins)
    #Offset the verts so we're in between the origin and maxs
    maxs-=mins
    for i in size(verts,1)
        verts[i]-=mins
    end
    
    #Find the longest dimensions
    gridSize, longestDim, ds = voxelMeshSize(maxs,resolution)
    origin = (1,1,1)
    xMax,yMax,zMax = gridSize
    grid = zeros(Bool,gridSize[1],gridSize[2],gridSize[3])

    #X-dir slices of array
    lines = getLines(ds,verts,tris)

    @assert(length(lines)%2 ==0)
    for i in 1:Int(length(lines)/2)
        p1 = lines[i*2-1]
        p2 = lines[i*2]

        x0 = clamp(round(Int64,p1[1]/ds),1,xMax)
        y0 = clamp(round(Int64,p1[2]/ds),1,yMax)
        z0 = clamp(round(Int64,p1[3]/ds),1,zMax)

        x1 = clamp(round(Int64,p2[1]/ds),1,xMax)
        y1 = clamp(round(Int64,p2[2]/ds),1,yMax)
        z1 = clamp(round(Int64,p2[3]/ds),1,zMax)
        
        @assert(x0==x1)
        if y0==y1 && z0 == z1
            grid[x0,y0,z0] = true
        else
            line(bresenhamf!,y0,z0,y1,z1,x0,grid)
        end
    end
    #Scanline fill 
    
    for i in 1:gridSize[1]
        for k in 1:gridSize[3]
            grid[i,:,k].=scanline(grid[i,:,k])
            
        end
    end
    return grid
end

mesh = load("C:/College/Research/AnisoEFITJulia/MeshFiles/IIWMiniASCII.stl");
#voxelizeMesh(mesh,resolution=100)
voxelizeMesh(mesh,resolution=100)
times = Vector{Float64}()
range = 50:100:5000
for i in range
    time = @elapsed res = voxelizeMesh(mesh, resolution=i)
    append!(times,time)
end
lines(range,times)
#GLMakie.volume(voxelizeMesh(mesh,resolution=100),isorange = 1e-10, isovalue=1, algorithm=:iso)