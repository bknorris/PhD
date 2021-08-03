function [ Data, Config ] = VPro_coordinateTransform( Data, Config, direction )
disp('Transforming velocity coordinates')
if nargin ~= 3
    disp( 'vectrinoIItransform requires three input arguments, the Config and Data structures as well as a transform direction' )
    return
end

direction = lower( direction );

if isfield( Config, 'nCells' )
    nCells = Config.nCells;
    nCalibratedCells = size( Config.ProbeCalibration_calibrationMatrix, 1 );
    originalPrefix = '';
elseif isfield( Config, 'Original_nCells' )
    nCells = Config.Original_nCells;
    nCalibratedCells = size( Config.ProbeCalibration_calibrationMatrix, 1 );
    originalPrefix = 'Original_';
end

switch direction
    case { 'xb' }
        if Config.( [ originalPrefix 'coordSystem' ] ) == 2
            disp( 'Already in beam coordinates' )
            return
        else
            Data.Profiles_VelBeam1 = NaN * zeros( size( Data.Profiles_VelX ) );
            Data.Profiles_VelBeam2 = NaN * zeros( size( Data.Profiles_VelX ) );
            Data.Profiles_VelBeam3 = NaN * zeros( size( Data.Profiles_VelX ) );
            Data.Profiles_VelBeam4 = NaN * zeros( size( Data.Profiles_VelX ) );
            for cell = 1:nCells
                if cell <= nCalibratedCells
                    T = reshape( Config.( [ originalPrefix 'ProbeCalibration_calibrationMatrix' ] )( cell, : ), 4, 4 )';
                else
                    T = eye( 4 );
                end
                V = double( [ Data.Profiles_VelX( :, cell )';
                    Data.Profiles_VelY( :, cell )';
                    Data.Profiles_VelZ1( :, cell )';
                    Data.Profiles_VelZ2( :, cell )' ] );
                B = inv( T ) * V;
                Data.Profiles_VelBeam1( :, cell ) = B( 1, : )';
                Data.Profiles_VelBeam2( :, cell ) = B( 2, : )';
                Data.Profiles_VelBeam3( :, cell ) = B( 3, : )';
                Data.Profiles_VelBeam4( :, cell ) = B( 4, : )';
            end
            Config.( [ originalPrefix 'coordSystem' ] ) = 2;
            Data = rmfield( Data, { 'Profiles_VelX', 'Profiles_VelY', 'Profiles_VelZ1', 'Profiles_VelZ2' } );
        end
        disp('Transformed XYZ to BEAM')
    case { 'bx' }
        if Config.( [ originalPrefix 'coordSystem' ] ) == 1
            disp( 'Already in XYZ coordinates.' )
            return
        else
            Data.Profiles_VelX = NaN * zeros( size( Data.Profiles_VelBeam1 ) );
            Data.Profiles_VelY = NaN * zeros( size( Data.Profiles_VelBeam1 ) );
            Data.Profiles_VelZ1 = NaN * zeros( size( Data.Profiles_VelBeam1 ) );
            Data.Profiles_VelZ2 = NaN * zeros( size( Data.Profiles_VelBeam1 ) );
            for cell = 1:nCells
                if cell <= nCalibratedCells
                    T = reshape( Config.( [ originalPrefix 'ProbeCalibration_calibrationMatrix' ] )( cell, : ), 4, 4 )';
                else
                    T = eye( 4 );
                end
                V = double( [ Data.Profiles_VelBeam1( :, cell )';
                    Data.Profiles_VelBeam2( :, cell )';
                    Data.Profiles_VelBeam3( :, cell )';
                    Data.Profiles_VelBeam4( :, cell )' ] );
                V = T * V;
                Data.Profiles_VelX( :, cell ) = V( 1, : )';
                Data.Profiles_VelY( :, cell ) = V( 2, : )';
                Data.Profiles_VelZ1( :, cell ) = V( 3, : )';
                Data.Profiles_VelZ2( :, cell ) = V( 4, : )';
            end
            Config.( [ originalPrefix 'coordSystem' ] ) = 1;
            Data = rmfield( Data, { 'Profiles_VelBeam1', 'Profiles_VelBeam2', 'Profiles_VelBeam3', 'Profiles_VelBeam4' } );
        end
        disp('Transformed BEAM to XYZ')
end