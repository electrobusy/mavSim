# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 16:09:01 2020

@author: Rohan
Supervisors: Guido de Croon & Christophe De Wagter 
"""

import numpy as np 

class flightgogglesv1():
    
    def __init__(self):
        
        # Simulation parameters: 
        self.dt_secs = 1/960
        
        # Convertion factors:
        self.deg2rad = np.pi/180

        # Vehicle properties:
        self.m = 1.0 # [kg]
        self.Inertia = np.array([[0.0049, 0, 0],[0, 0.0049, 0],[0, 0, 0.0069]]) # [kg m^2]
        self.momentArm = 0.08 # [m]
        self.motorTimeConstant = 0.02 # [sec]
        self.maxPropSpeed = 2200 # [rad/s] -- this information is per motor

        self.thrustCoeff = 1.91e-6 # [N/(rad/s)^2] -- C_Thrust
        self.torqueCoeff = 2.6e-7 # [N/(rad/s)^2]  -- C_torque
        self.dragCoeff = 0.1 #  [N/(m/s)^2]

        # Enviroment parameters:
        self.grav = 9.81 # [m/s^2]

        # Flags
        self.useAutoThrottle = False    
        self.includeDrag = True

        # Process Noise - Autocorrelation Value
        self.angAccelProcessNoiseAutoCorrelation = 0.00025 # [rad^2/s^2]
        self.linAccelProcessNoiseAutoCorrelation = 0.0005 # [m^2/s^3]   

        # Sensor Noise - Autocorrelation Value: 
        self.accMeasNoiseVariance = 0.005 # [m^2/s^4]
        self.gyroMeasNoiseVariance = 0.003 # [rad^2/s^2]
        
        # Filter the states - Low pass filter:
        self.gainP = 35530.5758439217
        self.gainQ = 266.572976289502

        # Inner loop - PID parameters: 
        self.propGain = np.array([9, 9, 9])
        self.intGain = np.array([3, 3, 3])
        self.derGain = np.array([0.3, 0.3, 0.3])
        self.intState = np.array([0, 0, 0])
        self.intBound = np.array([1000, 1000, 1000])
        
        # Bounding parameters: 
        self.CTRL_MAX_SPEED = 2.5
        self.CTRL_MAX_PITCH = 20*self.deg2rad
        self.CTRL_MAX_ROLL = 17*self.deg2rad
        self.CTRL_MAX_R = 20*self.deg2rad
        
        # -- States of the quadrotor (that and will vary throughout time and will be plotted):
        # low pass filter:
        self.filterState = np.zeros(3)
        self.filterStateDer = np.zeros(3)
        
        # low lever controller:
        self.intState = np.zeros(3)
        
        # Motor speeds: 
        self.propSpeedCommand = np.zeros(4)
        self.propSpeed = np.zeros(4)
        
        # Position/Velocity/Acceleration & Thrust & Drag
        self.position = np.zeros(3)
        self.velocity = np.zeros(3)
        self.specificForceBodyFrame = np.zeros(3)
        self.specificForce = np.zeros(3)
        
        self.specificThrust = 0
        
        # Attitude/AngVelocity/AngAcceleration
        self.attitude = np.zeros(4)
        self.attitudeDer = np.zeros(4)
        self.angVelocity = np.zeros(3)
        self.angAccel = np.zeros(3)
        
    # --- Functions: 
    def lpfproceedState(self, input_value):
        
        dt_secs = self.dt_secs
        
        # -- Determinant:
        det = self.gainP*dt_secs**2 + self.gainQ*dt_secs + 1

        # -- New filtered state derivative:
        stateDer = (self.filterStateDer + self.gainP*dt_secs*input_value)/det - dt_secs*self.gainP*self.filterState/det

        # -- New state:
        self.filterState = dt_secs*(self.filterStateDer + self.gainP*dt_secs*input_value)/det + ((dt_secs*self.gainQ + 1)*self.filterState)/det

        # -- Attribute the new filtered state derivative (for next iteration)
        self.filterStateDer = stateDer
        
        return
    
    def controlUpdate(self, command, curval, curder):
        
        dt_secs = self.dt_secs
        
        # - State deviation - error:
        stateDev = command - curval

        # - Integration of the state deviation - integral of the error
        self.intState = self.intState + dt_secs*stateDev
        
        # bound the integrative term 
        aux_conc1 = np.stack((-self.intBound, self.intState)) 
        aux_int = np.amax(aux_conc1, axis=0)
        
        aux_conc2 = np.stack((aux_int, self.intBound))
        self.intState = np.amin(aux_conc2, axis=0)

        # - PID control law:
        out = np.multiply(self.propGain,stateDev) + np.multiply(self.intGain, self.intState) + np.multiply(self.derGain, -curder)
        
        return out

    def computeMotorSpeedCommand(self, angAccCommand, thrust_z):

        # -- Compute momentThrust:
        momentThrust = np.concatenate((np.multiply(np.diag(self.Inertia),angAccCommand),thrust_z), axis=None)

        # -- Compute motorSpeedSquared:
        aux = 1/(4*self.momentArm)
        aux_C_T = aux/self.thrustCoeff
        aux_C_tau = aux/self.torqueCoeff
        G = np.array([
            [aux_C_T, -aux_C_T, -aux_C_tau*self.momentArm, aux_C_T*self.momentArm],
            [aux_C_T, aux_C_T, aux_C_tau*self.momentArm, aux_C_T*self.momentArm],
            [-aux_C_T, aux_C_T, -aux_C_tau*self.momentArm, aux_C_T*self.momentArm],
            [-aux_C_T, -aux_C_T, aux_C_tau*self.momentArm, aux_C_T*self.momentArm]])
        
        motorSpeedsSquared = G.dot(momentThrust)

        # -- Prop speed commands (assuming that rotor speeds >= 0) 
        aux1_conc = np.stack((np.zeros(motorSpeedsSquared.shape[0]),motorSpeedsSquared))
        aux1_max = np.sqrt(np.amax(aux1_conc, axis=0))
        
        aux2_conc = np.stack((aux1_max,self.maxPropSpeed*np.ones(motorSpeedsSquared.shape[0])))
        self.propSpeedCommand = np.amin(aux2_conc, axis=0)
        
        return
        
    
    def proceedState(self):
        
        dt_secs = self.dt_secs

        # Explicit Euler integration:
        # -- Position:
        self.position = self.position + dt_secs*self.velocity
        # -- Velocity:
        self.velocity = self.velocity + dt_secs*self.specificForce
        self.velocity[2] = self.velocity[2] - dt_secs*self.grav # Gravity influence

        # -- Attitude derivative: 
        # - Rotational tranformation matrix
        S_Q = np.array([
            [self.attitude[3], -self.attitude[2], self.attitude[1]],
            [self.attitude[2], self.attitude[3], self.attitude[0]],
            [-self.attitude[1], self.attitude[0], self.attitude[3]],
            [-self.attitude[0], -self.attitude[1], -self.attitude[2]]
            ])
        # - Compute attitude derivative in the quaternions coordinate frame:
        self.attitudeDer = 0.5*S_Q.dot(self.angVelocity);

        # Explicit Euler integration:
        # -- attitude:
        self.attitude = self.attitude + dt_secs*self.attitudeDer

        # -- Quaternions normalization factor:
        attNorm = np.sqrt(np.sum(self.attitude**2))

        # -- Normalize the quaternions:
        if attNorm > 1:  
            self.attitude = self.attitude/attNorm
        
        # -- Process noise:
        # - Angular Acceleration:
        angAccelProcessNoise = np.sqrt(self.angAccelProcessNoiseAutoCorrelation/dt_secs)*np.random.normal(np.zeros(3),np.ones(3))
        # -- Linear Acceleration:
        linAccelProcessNoise = np.sqrt(self.linAccelProcessNoiseAutoCorrelation/dt_secs)*np.random.normal(np.zeros(3),np.ones(3))

        #  Compute angular acceleration: -- use the Euler equation
        # self.angAccel = np.zeros(3)
        self.angAccel[0] = self.momentArm*self.thrustCoeff/self.Inertia[0,0]*(self.propSpeed[0]**2 + self.propSpeed[1]**2 - self.propSpeed[2]**2 - self.propSpeed[3]**2) + angAccelProcessNoise[0]
        self.angAccel[1] = self.momentArm*self.thrustCoeff/self.Inertia[1,1]*(-self.propSpeed[0]**2 + self.propSpeed[1]**2 + self.propSpeed[2]**2 - self.propSpeed[3]**2) + angAccelProcessNoise[1]
        self.angAccel[2] = self.torqueCoeff/self.Inertia[2,2]*(-self.propSpeed[0]**2 + self.propSpeed[1]**2 - self.propSpeed[2]**2 + self.propSpeed[3]**2) + angAccelProcessNoise[2]

        # Explicit Euler integration:
        # -- Angular velocity:
        self.angVelocity = self.angVelocity + dt_secs*self.angAccel
 
        # Propeller speed:
        aux_var = np.array([dt_secs, self.motorTimeConstant])
        aux_time = np.amin(aux_var)
        self.propSpeed = self.propSpeed + aux_time*((self.propSpeedCommand - self.propSpeed)/self.motorTimeConstant)
        
        # Vehicle specific thrust:
        self.specificThrust = self.thrustCoeff*(np.sum(self.propSpeed**2))/self.m

        # Vehicle total speed:
        speed = np.sqrt(np.sum(self.velocity**2))
        # print(speed)

        # [Added this myself - won't influence anything in the code]
        # dragForce = np.zeros(3)

        # Include the drag model: 
        if self.includeDrag and speed > 0:
            # Total drag:
            drag = self.dragCoeff*speed**2
            # Drag force:
            dragForce = drag*(self.velocity)/speed + linAccelProcessNoise

            # Express the drag in the body frame: 
            a = self.attitude[3]
            b = self.attitude[0]
            c = self.attitude[1]
            d = self.attitude[2]

            # Transformation matrix: 
            R = np.array([ 
                [a**2 + b**2 - c**2 - d**2, 2*b*c + 2*a*d, 2*b*d - 2*a*c],
                [2*b*c - 2*a*d, a**2 - b**2 + c**2 - d**2, 2*c*d + 2*a*b],
                [2*b*d + 2*a*c, 2*c*d - 2*a*b, a**2 - b**2 - c**2 + d**2],
                ])
            self.specificForceBodyFrame = R.dot(dragForce)
            self.specificForce = dragForce
        
        else:
            # Express the drag in the body frame: 
            a = self.attitude[3]
            b = self.attitude[0]
            c = self.attitude[1]
            d = self.attitude[2]

            # Transformation matrix: 
            R = np.array([ 
                [a**2 + b**2 - c**2 - d**2, 2*b*c + 2*a*d, 2*b*d - 2*a*c],
                [2*b*c - 2*a*d, a**2 - b**2 + c**2 - d**2, 2*c*d + 2*a*b],
                [2*b*d + 2*a*c, 2*c*d - 2*a*b, a**2 - b**2 - c**2 + d**2],
                ])
            
            self.specificForceBodyFrame = R.dot(linAccelProcessNoise)
            self.specificForce = linAccelProcessNoise

        self.specificForceBodyFrame[2] = self.specificForceBodyFrame[2] + self.specificThrust

        self.specificForce = self.specificForce + np.array([
            2*self.attitude[0]*self.attitude[2] + 2*self.attitude[1]*self.attitude[3],
            2*self.attitude[1]*self.attitude[2] - 2*self.attitude[0]*self.attitude[3],
            -self.attitude[0]**2 - self.attitude[1]**2 + self.attitude[2]**2 + self.attitude[3]**2,
            ])*self.specificThrust   

        return    

    def IMU_meas(self, angVel, accel):
        measured_angVel = angVel #  % + sqrt(self.gyroMeasNoiseVariance)*np.random.normal(np.zeros(3),np.ones(3))
        measured_accel = accel # % + sqrt(self.accMeasNoiseVariance)*np.random.normal(np.zeros(3),np.ones(3))

        measured_angVel_Cov = np.zeros(9)
        measured_acc_Cov = np.zeros(9)

        for i in range(0,9):
            if i == 1 or i == 5 or i == 9:
                measured_angVel_Cov[i] = self.gyroMeasNoiseVariance;
                measured_acc_Cov[i] = self.accMeasNoiseVariance;
            else:
                measured_angVel_Cov[i] = 0;
                measured_acc_Cov[i] = 0;
        
        return measured_angVel, measured_accel, measured_angVel_Cov, measured_acc_Cov
    
def Euler2Quat(phi, theta, psi):
    
    q = np.array([
         np.sin(phi/2)*np.cos(theta/2)*np.cos(psi/2) - np.cos(phi/2)*np.sin(theta/2)*np.sin(psi/2), # qx
         np.cos(phi/2)*np.sin(theta/2)*np.cos(psi/2) + np.sin(phi/2)*np.cos(theta/2)*np.sin(psi/2), # qy
         np.cos(phi/2)*np.cos(theta/2)*np.sin(psi/2) - np.sin(phi/2)*np.sin(theta/2)*np.cos(psi/2), # qz
         np.cos(phi/2)*np.cos(theta/2)*np.cos(psi/2) + np.sin(phi/2)*np.sin(theta/2)*np.sin(psi/2), # qw
         ])
    
    return q

def quat2Euler(q):
    #  -- Quaternions
    a = q[3]
    b = q[0]
    c = q[1]
    d = q[2]

    # -- Conversion from Quaternions to Euler:
    phi = np.arctan2(2*(a*b + c*d),1 - 2*(b**2 + c**2));
    theta = np.arcsin(2*(a*c - d*b));
    psi = np.arctan2(2*(a*d + b*c), 1 - 2*(c**2 + d**2));

    return np.array([phi, theta, psi])
        
if __name__ == '__main__':
        
    # -- Instantiate the FlightGoggles class 
    fg = flightgogglesv1()
    
    # -- simulation parameters: 
    # Simulation time 
    T = 4 # [sec]
    numPts = int(T/fg.dt_secs) # [-]
    # print(numPts)
    
    # -- Define initial states
    # initial attitude
    phi = 0 # [rad]
    theta = 0 # [rad]
    psi = 0 # [rad]
    
    fg.attitude = Euler2Quat(phi, theta, psi)
    
    # initial position: 
    fg.position = np.array([0, 0, 1.5]) # [m]
    
    # initial specific force (inertial frame): 
    fg.specificForce =  np.array([0, 0, -fg.grav])# [m s^-2]
    
    # initial propeller speed: 
    fg.propSpeed = np.sqrt(fg.m/4*fg.grav/fg.thrustCoeff)*np.ones(4) 
    
    # -- variables to plot later: 
    # Position/Velocity/Acceleration & Thrust & Drag
    position = list()
    velocity = list()
    specificForceBodyFrame = list()
    specificForce = list()
    
    specificThrust = list()
    
    # Attitude/AngVelocity/AngAcceleration
    attitude_q = list()
    attitude_E = list()
    angVelocity = list()
    angAccel = list()
    
    # Motor speeds: 
    propSpeedCommand = list()
    propSpeed = list()
    
    # -- inputs (angular velocity)
    # Note: Right now, only the inner loop rate controller is implemented. Thus,
    # only angular velocity reference commands are used.
    angVelCmd = np.zeros((numPts,3))
    thrust_z = fg.m*fg.grav
    
    # -- Simulation: 
    for i in range(0,numPts):
        # ----------- Things provided by Lockheed Martin --------------
        # 1 - Filter the measurements, in this case, the angular velocity:
        fg.lpfproceedState(fg.angVelocity)
        
        # 2 - Use the PID control law, in order to obtain the motor inputs:
        out = fg.controlUpdate(angVelCmd[i,:], fg.filterState, fg.filterStateDer)
        
        # 3 - Compute the motor speed commands, in order to obtain the desired response: 
        fg.computeMotorSpeedCommand(out, thrust_z)
       
        # 4 - Compute the state update for position, velocity, etc:
        fg.proceedState()
        
        # 5 - Obtain measurements from the IMU:
        fg.IMU_meas(fg.angVelocity, fg.specificForceBodyFrame)
        
        # ----------- Things provided by Lockheed Martin --------------
        
        # store variables for plotting: 
        position.append(fg.position)
        velocity.append(fg.velocity)
        specificForceBodyFrame.append(fg.specificForceBodyFrame)
        specificForce.append(fg.specificForce)
        
        specificThrust.append(fg.specificThrust)
        attitude_q.append(fg.attitude)
        angVelocity.append(fg.angVelocity)
        angAccel.append(fg.angAccel)
        
        EulerAng = quat2Euler(fg.attitude)
        attitude_E.append(EulerAng)
        
        propSpeedCommand.append(fg.propSpeedCommand)
        propSpeed.append(fg.propSpeed)
        
    # convert to numpy array for plotting
    position = np.array(position)
    velocity = np.array(velocity)
    specificForceBodyFrame = np.array(specificForceBodyFrame)
    specificForce = np.array(specificForce)
    
    specificThrust = np.array(specificThrust)
    attitude_q = np.array(attitude_q)
    angVelocity = np.array(angVelocity)
    angAccel = np.array(angAccel)
    
    attitude_E = np.array(attitude_E)
    
    propSpeedCommand = np.array(propSpeedCommand)
    propSpeed = np.array(propSpeed)
        
    # -- Plots:
    import matplotlib.pyplot as plt
    
    time = np.linspace(0,T,numPts)
    
    # position:
    fig, ax = plt.subplots(3,1)
    list_label_pos = [r'$x$ [m]', r'$y$ [m]', r'$z$ [m]']
    for i in range(0,3):
        ax[i].plot(time,position[:,i])
        ax[i].set_xlabel('t [sec]')
        ax[i].set_ylabel(list_label_pos[i])
        ax[i].grid()
    fig.tight_layout()
    fig.suptitle('Position')
    
    fig, ax = plt.subplots(3,1)
    list_label_vel = [r'$V_x$ [m/s]', r'$V_y$ [m/s]', r'$V_z$ [m/s]']
    for i in range(0,3):
        ax[i].plot(time,velocity[:,i])
        ax[i].set_xlabel('t [sec]')
        ax[i].set_ylabel(list_label_vel[i])
        ax[i].grid()
    fig.tight_layout()
    fig.suptitle('Velocity')
    
    fig, ax = plt.subplots(3,1)
    list_label_vel = [r'$A_x$ [m/s]', r'$A_y$ [m/s]', r'$A_z$ [m/s]']
    for i in range(0,3):
        ax[i].plot(time,specificForce[:,i])
        ax[i].set_xlabel('t [sec]')
        ax[i].set_ylabel(list_label_vel[i])
        ax[i].grid()
    fig.tight_layout()
    fig.suptitle('Specific Force')
    
    fig, ax = plt.subplots(3,1)
    list_label_vel = [r'$a_x$ [m/s]', r'$a_y$ [m/s]', r'$a_z$ [m/s]']
    for i in range(0,3):
        ax[i].plot(time,specificForceBodyFrame[:,i])
        ax[i].set_xlabel('t [sec]')
        ax[i].set_ylabel(list_label_vel[i])
        ax[i].grid()
    fig.tight_layout()
    fig.suptitle('Specific Force - Body Frame')
    
    fig, ax = plt.subplots(3,1)
    list_label_att = [r'$\phi$ [deg]', r'$\theta$ [deg]', r'$\psi$ [deg]']
    for i in range(0,3):
        ax[i].plot(time,attitude_E[:,i])
        ax[i].set_xlabel('t [sec]')
        ax[i].set_ylabel(list_label_att[i])
        ax[i].grid()
    fig.tight_layout()
    fig.suptitle('Attitude - Euler Angles')
    
    fig, ax = plt.subplots(2,2)
    list_label_att = [r'$\Omega_1$ [rpm]', r'$\Omega_2$ [rpm]', r'$\Omega_3$ [rpm]', r'$\Omega_4$ [rpm]']
    idx_row = [0, 0, 1, 1]
    idx_col = [0, 1, 0, 1]
    for i in range(0,4):
        ax[idx_row[i],idx_col[i]].plot(time,propSpeed[:,i],label='meas')
        ax[idx_row[i],idx_col[i]].plot(time,propSpeedCommand[:,i],label='cmd')
        ax[idx_row[i],idx_col[i]].set_xlabel('t [sec]')
        ax[idx_row[i],idx_col[i]].set_ylabel(list_label_att[i])
        ax[idx_row[i],idx_col[i]].grid()
    ax[0,0].legend()
    fig.tight_layout()
    fig.suptitle('Propeller Speeds')

        