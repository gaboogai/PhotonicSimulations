--length unit: um
--frequency: 1/wave length
--see the S4 doc
-- creates files : thickness_shape_lattice_radius.txt
-- files contain matrix with angle and wavelength as the rows and columns respectively

function round(num, numDecimalPlaces)
  local mult = 10^(numDecimalPlaces or 0)
  return math.floor(num * mult + 0.5) / mult
end

--------------------------------------------------------
--   1/ SET THE PARAMETERS
--------------------------------------------------------

-- wavelength range

lambda = 1.064

-- angle of incidence range

angle_start = 0
angle_stop = 1
n_angles = 1

-- thickness range

thickness = 0.02

-- pattern type

shape = "hex"

-- material

material = "SiN"

-- lattice parameter range

lattice_start = 1.21
lattice_stop = 1.23
n_lattices = 100

-- radius range

radius_start = 0.2
radius_stop = 0.6
n_radius = 100

--------------------------------------------------------
--   2/ REFLECTIVITY FUNCTION
--------------------------------------------------------

step_radius = (radius_stop - radius_start)/n_radius
step_lattice = (lattice_stop - lattice_start)/n_lattices
step_angle = (angle_stop - angle_start)/n_angles

tic = os.time()

for angle=angle_start,angle_stop,step_angle do

	filename = string.format("%s%s%s","cleanbroad",tostring(angle),".csv")
	f = io.open(filename, "w")

	f:write( ' ' .. ',' )

	for lattice=lattice_start,lattice_stop,step_lattice do

		f:write(lattice .. ',')	

	end

	f:write('\n')

	for radius=radius_start,radius_stop,step_radius do
	
		f:write(radius .. ',')			

		for lattice=lattice_start,lattice_stop,step_lattice do	
		
			if lattice<=2*radius then
		
				f:write( ' ' .. ',' )

			elseif lattice<0 then

				f:write( ' ' .. ',' )

			else

				S = S4.NewSimulation()	

				-- square lattice			
				-- S:SetLattice({lattice,0}, {0,lattice}) 

				-- hexagon lattice
				S:SetLattice({lattice,0}, {math.cos(math.pi/3)*lattice,math.sin(math.pi/3)*lattice})

				S:SetNumG(20)
			
				S:AddMaterial("Vacuum", {1,0})			
				S:AddMaterial("Si", {3.4422^2,0}) --(https://refractiveindex.info/?shelf=main&book=Si&page=Pierce 19/07/2018)
				S:AddMaterial("SiN", {1+3.0249*lambda^2/(lambda^2-0.1353406^2)+40314*lambda^2/(lambda^2-1239.842^2),0})	
				S:SetFrequency(1./lambda)
	
				S:AddLayer('AirAbove', 0 , 'Vacuum')
	
				S:AddLayer('Slab', thickness, material)

				S:AddLayer('AirBelow', 0, 'Vacuum')
	

				S:SetLayerPatternCircle('Slab',   -- layer to pattern
				'Vacuum', -- material inside circle
				{0,0},    -- center of circle
				radius)      -- radius of circle

				S:SetExcitationPlanewave(
				{angle,0}, -- incidence angles, phi in [0,180), theta in [0,360)
				{0,0}, -- s-polarization amplitude and phase (in degrees)
				{1,0}) -- p-polarization amplitude and phase


				forw_r, back_r, forw_i, back_i = S:GetPowerFlux('AirAbove')

				f:write(back_r .. ',')

			end
			
		end

		f:write('\n')

		print('completion: ' .. 100*(radius-radius_start+step_radius)/(radius_stop-radius_start) .. '%')

	end

end

toc = os.time()

print('Total time: ' .. (toc-tic)/60)
