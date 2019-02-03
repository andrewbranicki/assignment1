% PART 3: ENHANCEMENTS

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
    % Mean time between collisions is 0.2e-12
    MT_C = 0.2e-12;
    
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
    
    % THESE PARTICLES CANNOT BE IN THE ADDED RECTANGLES
    % Need to check their current starting location, and if it's inside
    % the rectangles then we need a new random position.
    
    % Check if particle is inside X boundary
    check_X_left = X_position > 0.8e-7;
    check_X_right = X_position < 1.2e-7;
    check_X = check_X_left & check_X_right;
    % Check which box it's in by checking Y coordinates
    check_top = Y_position > 0.6e-7;
    check_bottom = Y_position < 0.4e-7;
    % Find out the box
    box_top = check_top & check_X;
    box_bottom = check_bottom & check_X;
    IN_A_BOX = box_top | box_bottom;
    
    % Now we start randomizing the positions of particles inside box
    while(sum(IN_A_BOX) > 0)
        
        temp_x = rand(1,sum(IN_A_BOX));
        temp_y = rand(1,sum(IN_A_BOX));
        % Randomize position
        X_position(IN_A_BOX) = temp_x*region_x;
        Y_position(IN_A_BOX) = temp_y*region_y;
        %X_position(IN_A_BOX) = rand*region_x;
        %Y_position(IN_A_BOX) = rand*region_y;
        
        % Check again
        check_X_left = X_position > 0.8e-7;
        check_X_right = X_position < 1.2e-7;
        check_X = check_X_left & check_X_right;
        check_top = Y_position > 0.6e-7;
        check_bottom = Y_position < 0.4e-7;
        box_top = check_top & check_X;
        box_bottom = check_bottom & check_X;
        IN_A_BOX = box_top | box_bottom;
    end
    % Once we leave while loop, we have no particles in the boxes
    
    % Set up the angle for each particle in radians
    angle(1,:) = X(2,:)*2*pi;
    % This time we need to set up a Maxwell-Boltzmann distribution
    % for the velocity components and the average will be v_th
    sigma = sqrt(C.kb*T/C.m_n)/4;
    max_boltz_dist = makedist('Normal',v_th,sigma);
    velocity = random(max_boltz_dist,1,particles);
    figure(1)
    hist(velocity)
    title('Particle Velocity Histogram')
    X_velocity = v_change*velocity(1,:).*cos(angle(1,:));
    Y_velocity = v_change*velocity(1,:).*sin(angle(1,:));
    % Set up scattering percent
    PSCAT = 1 - exp(-v_change/MT_C);
    
    mfp_vec = zeros(1,particles);
    
    %% PAINTING THE CANVAS ------------------------------------------------
    
    % Go through loop for specified number of time steps
    for p = 1:timestep
        
        % Scattering electrons
        scatter = rand(1,particles);
        DID_I_SCATTER = scatter < PSCAT;
        % If PSCAT is greater than scatter, the particle scatters.
        % In this case, all the 1's in logical array DID_I_SCATTER
        % will be scattered and need to have recalculated stuff.
        
        % New random angle:
        angle(DID_I_SCATTER) = rand*2*pi;
        % New random velocity:
        velocity = random(max_boltz_dist,1,particles);
        X_velocity(DID_I_SCATTER) = v_change*velocity(DID_I_SCATTER) ...
            .*cos(angle(DID_I_SCATTER));
        Y_velocity(DID_I_SCATTER) = v_change*velocity(DID_I_SCATTER) ... 
            .*sin(angle(DID_I_SCATTER));
        
        % Keep track of how many particles scattered
        % to figure out mean free path
        mfp_vec(~DID_I_SCATTER) = mfp_vec(~DID_I_SCATTER)+step_size;
        % Anything that scattered will be set to 0
        mfp_vec(DID_I_SCATTER) = 0;
        
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
       
       % BOX REFLECTIONS -------------------------------------------------
       % BOTTOM BOX TEST
       % Check if particle is inside X boundary
        check_X_left = (X_position + X_velocity) > 0.8e-7;
        check_X_right = (X_position + X_velocity) < 1.2e-7;
        check_X = check_X_left & check_X_right;
        % Check Y boundary
        check_bottom = (Y_position + Y_velocity) < 0.4e-7;
        % Check bottom box
        box_bottom = check_bottom & check_X;
        X_velocity(box_bottom) = -1*X_velocity(box_bottom);
        % Check if we need to multiply Y component by -1 as well
        check_X_left = (X_position + X_velocity) > 0.8e-7 + step_size;
        check_X_right = (X_position + X_velocity) < 1.2e-7 - step_size;
        check_X = check_X_left & check_X_right;
        check_Y_bottombox_top = Y_position < 0.4e-7 - step_size;
        Y_BOX_BOTTOM = check_X & check_Y_bottombox_top;
        Y_velocity(Y_BOX_BOTTOM) = -1*Y_velocity(Y_BOX_BOTTOM);
        
        % TOP BOX TEST
        % Check if particle is inside X boundary
        check_X_left = (X_position + X_velocity) > 0.8e-7;
        check_X_right = (X_position + X_velocity) < 1.2e-7;
        check_X = check_X_left & check_X_right;
        % Check Y boundary
        check_top = (Y_position + Y_velocity) > 0.6e-7;
        % Check bottom box
        box_top = check_top & check_X;
        X_velocity(box_top) = -1*X_velocity(box_top);
        % Check if we need to multiply Y component by -1 as well
        check_X_left = (X_position + X_velocity) > 0.8e-7 + step_size;
        check_X_right = (X_position + X_velocity) < 1.2e-7 - step_size;
        check_X = check_X_left & check_X_right;
        check_Y_topbox_bottom = Y_position > 0.6e-7 + step_size;
        Y_BOX_TOP = check_X & check_Y_topbox_bottom;
        Y_velocity(Y_BOX_TOP) = -1*Y_velocity(Y_BOX_TOP);
        
       % Since we are using logical array, any array position with a '0'
       % means that the position did not go over the canvas boundary.
       % Therefore, we can just simply add the position and velocity
       % for these particles and not worry about the boundary stuff.
       
       % the cool particles are the ones that aren't in x_left or x_right
       % or in the boxes
       the_cool_particles = ~(X_left | X_right | box_top | box_bottom);
       X_position(the_cool_particles) = X_position(the_cool_particles) ...
           + X_velocity(the_cool_particles);
       
       % Since we're doing a whole bunch of time steps, we're gonna need
       % to save the current timestep results for future plotting.
       X_STORAGE(p,:) = X_position(1,:);
       Y_STORAGE(p,:) = Y_position(1,:);
       

       % Finally, calculate temperature
       avg_velocity = sum(velocity)/particles;
       current_temp(1,p) =  avg_velocity^2*C.m_n/(2*C.kb);
       % Mean free path calculation
       MFP = sum(mfp_vec)/particles;
       avg_tbc = MFP/avg_velocity;
       
    end
    
    %% SHOWING THE ART ----------------------------------------------------
    % Everything is calculated! Now let's start plotting...
    
    % Electron Density Map
    figure(2)
    electron_density = [X_position',Y_position'];
    hist3(electron_density,'CdataMode','auto')
    xlabel('X (m)')
    ylabel('Y (m)')
    title('Electron Density Map')
    colorbar
    view(2)
    
    % Temperature Map
    % Set up the X and Y linspace for the meshgrid
    X_canvas = linspace(0,region_x,10);
    Y_canvas = linspace(0,region_y,10);

    % Group x and y position values into bins defined in x and y canvas
    x_bin = discretize(X_position,X_canvas);
    y_bin = discretize(Y_position,Y_canvas);

    % Set up the temperature bins 
    TEMP_BINS = zeros(10,10);
    % Start going through the x and y values
    for x = 1:10
        for y=1:10
            % Define the x and y bins as logic values to easily assign
            % every value in the x and y velocities to a proper bin
            logicX = x_bin == x;
            logicY = y_bin == y;
            logic = logicX & logicY;
            % Whichever values in X and Y velocity fit into the correct bin
            % Will be added to be summed
            sum_x = sum(X_velocity(logic)) / v_change;
            sum_y = sum(Y_velocity(logic)) / v_change;
            
            velocity_avg = sqrt((sum_x)^2 + (sum_y)^2);
            % Use average velocity to calculate temperature
            % and put it into the proper bin
            TEMP_BINS(x,y) = velocity_avg^2*C.m_n/(2*C.kb);
        end
    end
    
    % Plot the bins as a temperature map
    figure(3)
    surf(TEMP_BINS)
    title('Electron Temperature Map')
    colorbar
    
    % Box Points for Plotting
    box_pic_x = [0.8e-7 0.8e-7  1.2e-7 1.2e-7];
    box_top_y = [1e-7 0.6e-7 0.6e-7 1e-7];
    box_bottom_y = [0 0.4e-7 0.4e-7 0];
    for row = 1:timestep
        figure(4)
        subplot(2,1,1);
        plot(X_STORAGE(row,:),Y_STORAGE(row,:),'o')
        xlim([0 region_x])
        ylim([0 region_y])
        xlabel('X (m)')
        ylabel('Y (m)')
        title('Electrons Flying Around Doing Cooler Things While Having Boxes in the Way')
        hold on
        plot(box_pic_x,box_top_y,'k')
        plot(box_pic_x,box_bottom_y,'k')
        subplot(2,1,2);
        plot(row,current_temp(row),'.k');
        title('The Current Temperature Caused by Those Electrons')
        ylabel('Temperature (K)')
        xlabel('Time-step')
        %  mean time between collisions is ?mn = 0.2ps
        legend(['Current Temperature:' num2str(current_temp(row))], ...
            ['Avg Time Between Collisions:' num2str(avg_tbc)], ...
            ['Mean Free Path:' num2str(MFP)])
        hold on
        xlim([0 timestep])
        ylim([250 350])
        pause(0.00000001)
    end