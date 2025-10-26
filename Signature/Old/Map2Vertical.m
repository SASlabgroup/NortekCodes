function [Data,zDist] = Map2Vertical(DataOrg, ConfigOrg)

Data = DataOrg;

plotData = false;

Row = ConfigOrg.burst_beam2xyz_rows;
Col = ConfigOrg.burst_beam2xyz_columns;
T_beam2xyz = transpose(reshape(ConfigOrg.burst_beam2xyz(1,1:Row*Col),Row,Col));

NB = ConfigOrg.burst_beam2xyz_columns;
NC = ConfigOrg.burst_nCells;
CS = ConfigOrg.burst_cellSize;
BD = ConfigOrg.burst_blankingDistance;
Nsamp=length(DataOrg.Burst_EnsembleCount);

if (findstr(ConfigOrg.sensorBoard,'AHRS') > 0 ),
   bAHRS = true;
else 
   bAHRS = false;
end

mappedCell = zeros(NC,NB);
zDist      = zeros(NC,NB);
fWeight1   = zeros(NC,NB);
fWeight2   = zeros(NC,NB);
nCell1     = zeros(NC,NB);
nCell2     = zeros(NC,NB);

% Save original data
Data.Burst_VelBeam1org = DataOrg.Burst_VelBeam1;
Data.Burst_VelBeam2org = DataOrg.Burst_VelBeam2;
Data.Burst_VelBeam3org = DataOrg.Burst_VelBeam3;
if NB == 4,
   Data.Burst_VelBeam4org = DataOrg.Burst_VelBeam4;
end


function i = c2m(j)
  i = j + 1;
end

function [valid] = validNC(nc, NC)
   valid = (nc >= 0) && (nc < NC);
end

function [v] = getBeam(n,c,b,type)
   v = DataOrg.(['Burst_' type 'Beam' char(48+b)])(n,c);
end

function setBeamVel(n,c,b,type,v)
   Data.(['Burst_' type 'Beam' char(48+b)])(n,c) = v;
end

function [beam2xyz3x3] = convBeam2XYZ(beam2xyz)
   beam2xyz3x3 = zeros(3,4);
   beam2xyz3x3(1,:) = beam2xyz(1,:); 
   beam2xyz3x3(2,:) = beam2xyz(2,:); 
   beam2xyz3x3(3,:) = (beam2xyz(3,:)+beam2xyz(4,:))/2;
end

if ConfigOrg.fwVersionDoppler <= 2163 && NB == 4,
    % Fix bug in transformation matrix in version 2163
    T_beam2xyz(2,:) = -T_beam2xyz(2,:);
end


% Change transformation matrix in the 4x4 case
if NB == 4
   M = pinv(convBeam2XYZ(T_beam2xyz));
else
   M=inv(T_beam2xyz);  % Desired cell position
end


% Precalculate along beam cell size and blanking distance
CSslant = CS/M(c2m(0),c2m(2)); % Correct for slant angle
BDslant = BD/M(c2m(0),c2m(2)); % Correct for slant angle
slantFactor = M(c2m(0),c2m(2));

% Precalculate zDist
for i=0:NC-1,
   for j=0:NB-1,
      zDir  = M(c2m(j),c2m(2));
      % Very general calculations, more simply zDist = BD+CS*(i+1) could be
      % used
      zDist(c2m(i),c2m(j)) = zDir*(BDslant + CSslant + i*CSslant); 
   end
end
 
% Change transformation matrix in the 4x4 case
if NB == 4
   [T_beam2xyz] = convBeam2XYZ(T_beam2xyz);
end

T_beam2xyzCurrent = T_beam2xyz;

if Nsamp > 1,
   h = waitbar(0,'Processing...');
end

for k=0:Nsamp-1,
   % Calculate vertical projection for the various beams
   if bAHRS
      xyz2enu = transpose(reshape(DataOrg.Burst_AHRSRotationMatrix(c2m(k),:),3,3));
   else
      xyz2enu = CalcHPRMatrix(DataOrg.Burst_Heading(c2m(k)), DataOrg.Burst_Pitch(c2m(k)), DataOrg.Burst_Roll(c2m(k)));
      if (bitand(bitshift(uint32(DataOrg.Burst_Status(c2m(k))), -25),7) == 5)
         T_beam2xyzCurrent(1,:) = T_beam2xyz(1,:);
         T_beam2xyzCurrent(2:end,:) = -T_beam2xyz(2:end,:);
      else
         T_beam2xyzCurrent = T_beam2xyz;
      end
   end 

   
   
   if NB == 4 
      M = pinv(xyz2enu*T_beam2xyzCurrent);
   else
      M = inv(xyz2enu*T_beam2xyzCurrent);
   end
   
   for i=0:NC-1,
      for j=0:NB-1,
         % Determine the matching cell for the tilted beam
         zDirM  = abs(M(c2m(j),c2m(2)));
         
         % The projected vertical distance is given by: zDirM*(BDslant + (i+1)*CSslant)
         % Solved for the desired vertical distance zDist this gives the
         % mapped cell position
         mappedCell(c2m(i),c2m(j)) = -1 + ((zDist(c2m(i))/zDirM) - BDslant)/CSslant;
         fMapdCell = mappedCell(c2m(i),c2m(j));
         if fMapdCell > -0.1 && fMapdCell < 0.0
            fMapdCell = 0.0;
         end
                  
         if fMapdCell >= 0.0 && fMapdCell < NC
            % Use neighbours to weight linearly
            fWeight1(c2m(i),c2m(j)) = ceil(fMapdCell) - fMapdCell;
            fWeight2(c2m(i),c2m(j)) = fMapdCell   - floor(fMapdCell);
            nCell1(c2m(i),c2m(j))   = floor(fMapdCell);
            nCell2(c2m(i),c2m(j))   = ceil(fMapdCell);

            if  fWeight1(c2m(i),c2m(j)) == 0.0 && fWeight2(c2m(i),c2m(j)) == 0.0,
                fWeight1(c2m(i),c2m(j)) = 1.0;
                nCell1(c2m(i),c2m(j))   = fMapdCell;
            end        

            nC1 = round(nCell1(c2m(i),c2m(j)));
            nC2 = round(nCell2(c2m(i),c2m(j)));

            if (validNC(nC1, NC) &&  validNC(nC2, NC))
               % Calculate interpolated velocity               
               vel = fWeight1(c2m(i),c2m(j))*getBeam(c2m(k),c2m(nC1),c2m(j),'Vel')+ ... 
                     fWeight2(c2m(i),c2m(j))*getBeam(c2m(k),c2m(nC2),c2m(j),'Vel');
               % Store the bin mapped velocity for this bin   
               setBeamVel(c2m(k),c2m(i),c2m(j),'Vel', vel);

               % Interpolate amplitude
               amp = fWeight1(c2m(i),c2m(j))*getBeam(c2m(k),c2m(nC1),c2m(j),'Amp')+ ... 
                     fWeight2(c2m(i),c2m(j))*getBeam(c2m(k),c2m(nC2),c2m(j),'Amp');
               setBeamVel(c2m(k),c2m(i),c2m(j),'Amp', amp);
               
               % Interpolate correlation
               cor = fWeight1(c2m(i),c2m(j))*getBeam(c2m(k),c2m(nC1),c2m(j),'Cor')+ ... 
                     fWeight2(c2m(i),c2m(j))*getBeam(c2m(k),c2m(nC2),c2m(j),'Cor');
               setBeamVel(c2m(k),c2m(i),c2m(j),'Cor', cor);                        
            else
               % Set velocity to invalid status and set correlation to 0% to
               % discard this velocity from averaging               
               setBeamVel(c2m(k),c2m(i),c2m(j),'Vel', -32.768);
               setBeamVel(c2m(k),c2m(i),c2m(j),'Cor', 0);               
            end
         else 
               % Set velocity to invalid status and set correlation to 0% to
               % discard this velocity from averaging
               setBeamVel(c2m(k),c2m(i),c2m(j),'Vel', -32.768);
               setBeamVel(c2m(k),c2m(i),c2m(j),'Cor', 0);               
         end
      end % j
   end % i
   if Nsamp > 1,
      waitbar(k/Nsamp,h)
   end
end % k

if Nsamp > 1,
   close(h)
end

if plotData
%    disp(M)
   nb=1;
   figure(nb)
   subplot(2,1,1)
   plot([nCell1(:,nb) nCell2(:,nb) mappedCell(:,nb) ],'.-'),grid
   legend('nCell1','nCell2','mappedCell','Location','SouthEast')
   subplot(2,1,2)
   plot([fWeight1(:,nb) fWeight2(:,nb)],'.-'),grid
   legend('fWeight1','fWeight2','Location','WestOutside')
   nb=2;
   figure(nb)
   subplot(2,1,1)
   plot([nCell1(:,nb) nCell2(:,nb) mappedCell(:,nb) ],'.-'),grid
   legend('nCell1','nCell2','mappedCell','Location','SouthEast')
   subplot(2,1,2)
   plot([fWeight1(:,nb) fWeight2(:,nb)],'.-'),grid
   legend('fWeight1','fWeight2','Location','WestOutside')
   nb=3;
   figure(nb)
   subplot(2,1,1)
   plot([nCell1(:,nb) nCell2(:,nb) mappedCell(:,nb) ],'.-'),grid
   legend('nCell1','nCell2','mappedCell','Location','SouthEast')
   subplot(2,1,2)
   plot([fWeight1(:,nb) fWeight2(:,nb)],'.-'),grid
   legend('fWeight1','fWeight2','Location','WestOutside')

   if NB == 4
      nb=4;
      figure(nb)
      subplot(2,1,1)
      plot([nCell1(:,nb) nCell2(:,nb) mappedCell(:,nb) ],'.-'),grid
      legend('nCell1','nCell2','mappedCell','Location','SouthEast')
      subplot(2,1,2)
      plot([fWeight1(:,nb) fWeight2(:,nb)],'.-'),grid
      legend('fWeight1','fWeight2','Location','WestOutside')
   end
end
%

end
