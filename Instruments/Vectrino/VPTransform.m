function [tdat1,tdat2,tdat3,tdat4] = VPTransform(dat1,dat2,dat3,dat4,TMatrix,direction)
%Input: dat1-dat4 - velocites either in BEAM or XYZ (where Z is Z1 and Z2
%       for a total of four fields.
%       TMatrix - transformation matrix for the specific Vectrino (found in
%       Config struct under the subfield ProbeCalibration_calibrationMatrix
%       direction - direction of transformation, either 'bx' for 'beam to
%       xyz' or 'xb' for 'xyz to beam'
%
%Output: tdat1-tdat4 - transformed velocities into the coordinate system of
%        choosing.
%

disp('Transforming velocity coordinates')

direction = lower( direction );
nCells = 35;
nCalibratedCells = size(TMatrix,1);

switch direction
    case { 'xb' }
        tdat1 = NaN*zeros(size(dat1));
        tdat2 = NaN*zeros(size(dat1));
        tdat3 = NaN*zeros(size(dat1));
        tdat4 = NaN*zeros(size(dat1));
        for cell = 1:nCells
            if cell <= nCalibratedCells
                T = reshape(TMatrix(cell,:),4,4)';
            else
                T = eye(4);
            end
            V = double([dat1(:,cell)';...
                dat2(:,cell)';...
                dat3(:,cell)';...
                dat4(:,cell)']);
            B = inv(T)*V;
            tdat1(:,cell) = B(1,:)';
            tdat2(:,cell) = B(2,:)';
            tdat3(:,cell) = B(3,:)';
            tdat4(:,cell) = B(4,:)';
        end
        disp('Transformed XYZ to BEAM')
    case { 'bx' }
        tdat1 = NaN*zeros(size(dat1));
        tdat2 = NaN*zeros(size(dat1));
        tdat3 = NaN*zeros(size(dat1));
        tdat4 = NaN*zeros(size(dat1));
        for cell = 1:nCells
            if cell <= nCalibratedCells
                T = reshape(TMatrix(cell,:),4,4)';
            else
                T = eye(4);
            end
            V = double([dat1(:,cell)';...
                dat2(:,cell)';...
                dat3(:,cell)';...
                dat4(:,cell)']);
            V = T*V;
            tdat1(:,cell) = V(1,:)';
            tdat2(:,cell) = V(2,:)';
            tdat3(:,cell) = V(3,:)';
            tdat4(:,cell) = V(4,:)';
        end
        disp('Transformed BEAM to XYZ')
end