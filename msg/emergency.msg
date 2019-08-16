uint64 timestamp			# time since system start (microseconds)

uint16 emergency_type
float32[3] n_axis
bool n_axis_valid			# Set to true if n_axis is valid
float32[3] axis_sp
float32[3] other_att_control
float32[3] other_att_control2
float32[3] att_control
float32[3] eul
float32[3] eul_sp
float32[3] rates
float32[3] rates_sp

uint16 NO_BROKEN_ROTORS = 0
uint16 ONE_BROKEN_ROTOR = 1
uint16 TWO_BROKEN_ROTORS = 2
uint16 THREE_BROKEN_ROTORS = 3

# TOPICS emergency emergency_info emergency_info2 lqr_control pid_control
