PLANET_MASS 5.2915793E22
PLANET_RADIUS 600000
PLANET_SCALE_HEIGHT 5000
PLANET_P0 1
PLANET_ROTATION_PERIOD 21600
PLANET_SOI 84159286
LAUNCH_LATITUDE -0.001691999747072392
LAUNCH_LONGITUDE 0
LAUNCH_ALTITUDE 77.6
MAX_VELOCITY 10000
NAME Kasuha_Cargo_to_Orbit_IV
ADD_STAGE 532.488 32 0.2
ADD_ENGINE 175 388 390 54
ADD_ENGINE 120 220 800
ADD_ENGINE  60 220 800
ADD_ENGINE  60 220 800
ADD_STAGE 490.434 32 0.2
ADD_ENGINE 175 388 390 50
ADD_ENGINE  60 220 800 4
ADD_STAGE 448.380 32 0.2
ADD_ENGINE 175 388 390 46
ADD_ENGINE 240 220 800
ADD_STAGE 406.326 32 0.2
ADD_ENGINE 175 388 390 42
ADD_ENGINE 240 220 800
ADD_STAGE 364.272 32 0.2
ADD_ENGINE 175 388 390 38
ADD_ENGINE  240 220 800
ADD_STAGE 322.218 32 0.2
ADD_ENGINE 175 388 390 34
ADD_ENGINE 240 220 800
ADD_STAGE 280.164 32 0.2
ADD_ENGINE 175 388 390 30
ADD_ENGINE 240 220 800
ADD_STAGE 229.100 32 0.2
ADD_ENGINE 175 388 390 20
ADD_ENGINE  240 220 800
ADD_STAGE 178.036 32 0.2
ADD_ENGINE 175 388 390 10
ADD_ENGINE 240 220 800
ADD_STAGE 126.972 16 0.2
ADD_ENGINE 240 220 800
TARGET_PERIAPSIS 75000
ITERATIONS 50
SET_NODES 10 20
MESH_REFINEMENT manual
NLP_TOLERANCE 1.0e-5
COMPUTE
#TEST
POSTPROCESS
#GET_CONTROLS 0 600 50
#GET_CONTROLS 122
#GET_PITCH_THRUST 10
#GET_PITCH_THRUST 100
#GET_PITCH_THRUST 700
GET_FINAL_TIMES
GET_PITCH_THRUST 10
