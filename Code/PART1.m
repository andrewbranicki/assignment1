% PART 1 ELECTRON MODELLING

particles = 1000;

    global C X Y
    C.q_0 = 1.60217653e-19;             % electron charge
    C.hb = 1.054571596e-34;             % Dirac constant
    C.h = C.hb * 2 * pi;                % Planck constant
    C.m_0 = 9.10938215e-31;             % electron mass
    C.kb = 1.3806504e-23;               % Boltzmann constant
    C.eps_0 = 8.854187817e-12;          % vacuum permittivity
    C.mu_0 = 1.2566370614e-6;           % vacuum permeability
    C.c = 299792458;                    % speed of light
    C.g = 9.80665;                      % metres (32.1740 ft) per s²
    C.m_n = 0.26*C.m_0;                 % effective mass of electrons

    %% SETTING UP THE CANVAS ----------------------------------------------
    
    % The nominal size of the region is 200nm × 100nm
    region_x = 200e-9;
    region_y = 100e-9;
    % Pick time step size which should be smaller than 1/100 the region
    step_size = 1e-9;
    timestep = 1000;
    % Assume T = 300K
    T = 300;
    % Calculate thermal velocity
    v_th = sqrt(2*C.kb*T/C.m_n);
    % Change in velocity for each timestep
    v_change = step_size/v_th;
    
    %% RANDOM PARTICLE INITIALIZATION -------------------------------------
    
    % Creates 2 rows of random numbers for X coordinates of each particle.
    % The two rows are used for position and angle.
    % Y coordinates only requires position.
    X = rand(2,particles);
    Y = rand(1,particles);
    % Set up the X and Y position coordinates for each particle
    % Multiply by region restrictions to ensure it's within the region
    X_position(1,:) = X(1,:)*region_x;
    Y_position(1,:) = Y(1,:)*region_y;
    % Set up the angle for each particle in radians
    angle(1,:) = X(2,:)*2*pi;
    % Set up X and Y components of initial particle velocity using v_th
    % Shout out to grade 12 physics :^)
    X_velocity = v_th*v_change*cos(angle(1,:));
    Y_velocity = v_th*v_change*sin(angle(1,:));
    
    
    %% PAINTING THE CANVAS ------------------------------------------------
    
    % Go through loop for specified number of time steps
    for p = 1:timestep
        
        % Begin moving the electrons by using initial position
        % and adding the velocity to it
        % Use logical array to make stuff easier...
        
        % Keep all the particles within the canvas.
        % If their position + velocity goes over the right X boundary,
        % We will subtract region_x to make it jump to the opposite edge.
        X_right = X_position + X_velocity > region_x;
        X_position(X_right) = X_position(X_right) ... 
            + X_velocity(X_right) - region_x;
        % Do the same thing for the left X boundary (x = 0)
        X_left = X_position + X_velocity < 0;
        X_position(X_left) = X_position(X_left) ... 
            + X_velocity(X_left) + region_x;
        
       % Since we are using logical array, any array position with a '0'
       % means that the position did not go over the canvas boundary.
       % Therefore, we can just simply add the position and velocity
       % for these particles and not worry about the boundary stuff.
       
       % the cool particles are the ones that aren't in x_left or x_right
       the_cool_particles = ~(X_left | X_right);
       X_position(the_cool_particles) = X_position(the_cool_particles) ...
           + X_velocity(the_cool_particles);
       
       % Now we need to worry about the Y boundary...
       % For this, we can check if the Y coordinate is over the limit,
       % and if it is we will just make the velocity negative to turn
       % the particle around!
       % Check if the Y coordinate isn't cool:
       is_Y_not_cool = (Y_position + Y_velocity > region_y ... 
           | Y_position + Y_velocity < 0);
       % For the Y particles that aren't cool (1's in the logical array),
       % we will need to flip their velocities and turn them onto the path
       % of cool.
       Y_velocity(is_Y_not_cool) = -1*Y_velocity(is_Y_not_cool);
       % Now that the particles are all cool, we can update the positions
       Y_position(1,:) = Y_position(1,:) + Y_velocity(1,:);
       
       % Since we're doing a whole bunch of time steps, we're gonna need
       % to save the current timestep results for future plotting.
       X_STORAGE(p,:) = X_position(1,:);
       Y_STORAGE(p,:) = Y_position(1,:);
       
       % Finally, calculate temperature to make sure it stays at 300!
       SUM_X = sum( (X_velocity/v_change).^2);
       SUM_Y = sum( (Y_velocity/v_change).^2);
       temp = (SUM_X+SUM_Y)*C.m_n/(2*C.kb);
       current_temp(p) = temp/particles;
    end
    
    %% SHOWING THE ART ----------------------------------------------------
    % Everything is calculated! Now let's start plotting...
    figure(1)
    for paths = 1:particles
        subplot(3,1,1);
        plot(X_STORAGE(:,paths),Y_STORAGE(:,paths),'-')
        xlim([0 region_x])
        ylim([0 region_y])
        xlabel('X (m)')
        ylabel('Y (m)')
        title('Final Electron Paths')
        hold on
    end
    for row = 1:timestep
        subplot(3,1,2);
        plot(X_STORAGE(row,:),Y_STORAGE(row,:),'o')
        xlim([0 region_x])
        ylim([0 region_y])
        xlabel('X (m)')
        ylabel('Y (m)')
        title('Electrons Flying Around Doing Cool Things')
        hold on
        subplot(3,1,3);
        plot(row,current_temp(row),'.k');
        title('The Current Temperature Caused by Those Electrons')
        ylabel('Temperature (K)')
        xlabel('Time-step')
        %  mean time between collisions is ?mn = 0.2ps
        legend(['Current Temperature:' num2str(current_temp(row))], ...
            ['Thermal Velocity:' num2str(v_th)], ...
            ['Mean Free Path:' num2str(v_th*0.2e-12)])
        hold on
        xlim([0 timestep])
        ylim([T-1 T+1])
        pause(0.01)
    end
