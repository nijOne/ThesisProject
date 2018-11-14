classdef Allegretto
                                         
    properties % Variables acquired from PDF file

        image;
        eyeID;
        correctionVal;
        clinicalVal;
        targetVal;
        opticalZone;
        transitionZone;
        ablationZone;
        vertexDistance;
        K1;
        K1Axis;
        K2;
        K2Axis;
        DAx;
        pupilDiameter;
        maxDepth;
        centralDepth;
        cornealThickness;
        flapThickness;
        stroma;
        memo;

    end
    
    properties % Variables necessary for calculations
        
        x;
        y;
        X;
        Y;
        eps;
        border;
        K1Post;
        K2Post;
        edgeWidth;
        abZ;
        abZEdge;
        abTrans;
        zPre;
        zPost;
        R1;
        R2;
        R1Post;
        R2Post;
        rotAngle;
        grayPix;
        pixelTab;
        valueTab;
        scale;
        initImgSize;
        numOfColors;
        z;
        grid;
    end
    
    methods
        
        function obj = Allegretto(filename)
            
            savename = erase(filename, '.pdf'); 
            loadname = [replace(savename, 'ABLACJE', 'MAT_FILES') '.mat']
            
            if exist(loadname, 'file')
              load(loadname);
            end
            
                initialize(obj);
                readPDF(obj, filename);
                setProps(obj, 0.025, 5, 0.1);                
                obj = run(obj); 
                obj = calcABZ(obj);
                obj = calcEdge(obj);
                obj = calcTrans(obj);             
                save(savename);
  
        end
        
        function obj = initialize(obj)
                       
                obj.initImgSize = length(obj.image);
                obj.pixelTab = imread('pixTab.png');
                obj.grayPix = zeros(1,1,3,'uint8');
                obj.grayPix(1,1,1) = 160;
                obj.grayPix(1,1,2) = 160;
                obj.grayPix(1,1,3) = 160;
                obj.numOfColors = length(obj.pixelTab);
                obj.valueTab = round(linspace(obj.maxDepth, 0, obj.numOfColors),1);
                obj.valueTab(25) = NaN;
                obj.grid = zeros(1, obj.numOfColors);
                
        end
        
        function obj = calcABZ(obj)
            
                obj = zPreFunc(obj);
                obj = zPostFunc(obj);
                
                obj.abZ = (obj.zPre - obj.zPost)*1000;
                obj.abZ(obj.abZ<0) = 0;
                
        end
        
        function obj = setProps(obj, myEps, myBorder, myEdgeWidth)
            
                obj.eps = myEps;
                obj.border = myBorder;
                obj.edgeWidth = myEdgeWidth;
                
                obj.x = -obj.border:obj.eps:obj.border; 
                obj.y = -obj.border:obj.eps:obj.border; 
                [obj.X, obj.Y] = meshgrid(obj.x, obj.y);
                
        end
        
        function obj = run(obj)
           scaling(obj);
          
           sumGray = sum(obj.image(:) == obj.grayPix); 
           
           while ne(sumGray(1), 0)
               obj = clearImg(obj);
               sumGray = sum(obj.image(:) == obj.grayPix);
           end
           
           obj = processImage(obj);
           obj = calcKPost(obj);
           obj = calcR(obj);
                
           obj.R1 = obj.R1*1000; 
           obj.R2 = obj.R2*1000;
           obj.R1Post = obj.R1Post*1000;
           obj.R2Post = obj.R2Post*1000;
           
        end
        
        function DAx = distribution (obj, ax)

            DAx = obj.correctionVal(1) + cos((ax-obj.correctionVal(3))*pi/180)^2*obj.correctionVal(2);
            DAx = round(DAx*4)/4;
        end
        
        function obj = zPreFunc(obj)

            tempZ = zeros(length(obj.x),length(obj.y));
            
            for i = 1 : length(obj.x)
                for j = 1 : length(obj.y)
                    
                    if obj.x(i)^2 + obj.y(j)^2 <= (obj.ablationZone/2)^2
                        tempZ(i,j) = real(obj.R2.*sqrt(1 - (obj.x(i).^2 ./ obj.R1.^2 +obj. y(j).^2 ./ obj.R2.^2)));
                    else
                        tempZ(i,j) = NaN;
                    end
                    
                end
            end
            
            obj.zPre = tempZ;              
         end

        function obj = zPostFunc(obj)

            diff = obj.R2 - obj.R2Post;
            
            tempZ = zeros(length(obj.x),length(obj.y));
          
            R = obj.opticalZone;
        
            for i = 1 : length(obj.x)
                for j = 1 : length(obj.y)
                    
                    if obj.x(i)^2 + obj.y(j)^2 <= (R/2)^2
                        tempZ(i,j) = real(obj.R2Post.*sqrt(1 - (obj.x(i).^2 ./ obj.R1Post.^2 + obj.y(j).^2 ./ obj.R2Post.^2))) + diff - obj.centralDepth*0.001;
                    else
                        tempZ(i,j) = NaN;
                    end
                    
                end
            end
            
            obj.zPost = tempZ;            
        end
        
        function nearest = getNEdge(obj,Xp,Yp)
            
            Z = obj.abZEdge;
            maska = not(isnan(Z));
            mX = obj.X(maska);
            mY = obj.Y(maska);
            
            dist = sqrt((mX - Xp).^2 + (mY - Yp).^2);
            
            [~, id] = min(dist);
            
            xIdx = find (obj.x == mX(id));
            yIdx = find (obj.y == mY(id));
            
            nearest = obj.abZEdge(xIdx,yIdx);
           
        end
        
        function obj = calcEdge(obj)

            tempEdge = zeros(length(obj.x),length(obj.y));

            for i = 1 : length(obj.x)
                for j = 1 : length(obj.y)
                    
                    if obj.x(i)^2 + obj.y(j)^2 >= ((obj.opticalZone/2 - obj.edgeWidth))^2 && obj.x(i)^2 + obj.y(j)^2 <= (obj.opticalZone/2)^2
                        tempEdge(i,j) = obj.abZ(i,j);
                    else
                        tempEdge(i,j)  = NaN;
                    end                  
                end
            end
            
            obj.abZEdge = tempEdge;
          
        end
  
        function obj = calcTrans(obj)
            
            RR = obj.opticalZone + 2.*obj.transitionZone;
            
            tempTrans = zeros(length(obj.x),length(obj.y));
            
            for i = 1 : length(obj.x)
                
                for j = 1 : length(obj.y)
                    
                    if (obj.x(i)^2 + obj.y(j)^2 >= (obj.opticalZone/2)^2 && obj.x(i)^2 + obj.y(j)^2 <= ((RR)/2)^2)
                        if  not(isnan(getNEdge(obj,obj.x(i),obj.y(j)))) && getNEdge(obj,obj.x(i),obj.y(j)) >= 0.02
                        tempTrans(i,j) = getNEdge(obj,obj.x(i),obj.y(j)) * (1 + cos(pi/2*((sqrt(obj.x(i)^2 + obj.y(j)^2) - 0.5*obj.opticalZone)/(0.5*RR - 0.5*obj.opticalZone) + 1)));
                        obj.abZ(i,j) = tempTrans(i,j);
                        end
                    else
                        tempTrans(i,j) = NaN;
                    end
                    
                end 
            end 
            
            obj.abTrans = tempTrans;
            
        end
        
        function myEllis(obj)

                mZPre = min(obj.zPre(:));
                mZPost = min(obj.zPost(:));
                ZPre = obj.zPre;% - mZPre ;
                ZPost = obj.zPost;%  - mZPost;% + (obj.R2 - obj.R2Post);
                
                 %hold on
                %surf(obj.x,obj.y,ZPre, 'FaceColor', [1 0 0],'FaceAlpha', 0.5,'EdgeAlpha', 0.0);
                %surf(obj.x,obj.y,ZPost, 'FaceColor', [0 1 0],'FaceAlpha', 1,'EdgeAlpha', 0.0);
                mesh(obj.x,obj.y,obj.abZ);% 'FaceColor', [0 1 0],'FaceAlpha', 1,'EdgeAlpha', 0.1);
                xlabel('x[mm]','FontSize', 20);
                ylabel('y[mm]','FontSize', 20);
                zlabel('ablation depth[um]','FontSize', 20);
                 %axis equal
                
                 hold off
                      
        end
        
        function obj = fillProps(obj, text)

            c1 = @(str)textscan(str{1}, '%f%*s');
            c2 = @(str)textscan(str{1}, ' %*f D @ %f°');
            c3 = @(str)textscan(str{1}, ' %*s %f %*s');
            c4 = @(str)textscan(str{1}, ' %s');
            
            getD = @(c)c{1};    
            obj.eyeID = getD(c4(text(6)));
            obj.eyeID = obj.eyeID{1};
            obj.correctionVal = regexp(erase(erase(text{12}, '@'),'°'), 'D\s', 'split');
            obj.correctionVal = sscanf(sprintf('%s', obj.correctionVal{:}),'%f');
            obj.clinicalVal = regexp(erase(erase(erase(text{13}, '@'),'°'),'^'), 'D\s', 'split');
            obj.clinicalVal = sscanf(sprintf('%s', obj.clinicalVal{:}),'%f');        
            obj.targetVal = regexp(erase(erase(erase(text{14}, '@'),'°'),'^'), 'D\s', 'split');
            obj.targetVal = sscanf(sprintf('%s', obj.targetVal{:}),'%f');
            obj.opticalZone = getD(c1(text(18)));
            obj.transitionZone = getD(c1(text(19)));
            obj.ablationZone = getD(c1(text(20)));
            obj.vertexDistance = getD(c1(text(24)));
            obj.K1 = getD(c1(text(25)));           
            obj.K1Axis = getD(c2(text(25)));
            obj.K2 = getD(c1(text(26)));
            obj.K2Axis = getD(c2(text(26)));
            obj.pupilDiameter = getD(c1(text(29)));
            obj.maxDepth = getD(c1(text(39)));
            obj.centralDepth = getD(c1(text(40)));
            obj.cornealThickness = getD(c1(text(44)));
            obj.flapThickness = getD(c1(text(45)));
            obj.stroma = getD(c3(text(46)));
            
            if ne(string(text(48)), 'WaveLight AG')
                obj.memo = text(48);
            end
            
        end
        
        function obj = readPDF(obj, filename)

            javaaddpath('pdfbox-1.8.3.jar');
            import org.apache.pdfbox.cos.COSName;
            import org.apache.pdfbox.pdmodel.PDDocument;
            import org.apache.pdfbox.pdmodel.PDPage;
            import org.apache.pdfbox.pdmodel.PDPageTree;
            import org.apache.pdfbox.pdmodel.PDResources;
            import org.apache.pdfbox.text.PDFTextStripper;
            import org.apache.pdfbox.pdmodel.graphics.PDXObject;
            import org.apache.pdfbox.pdmodel.graphics.image.PDImageXObject;

            pdfName = java.lang.String(filename);
            pdfDoc = PDDocument;
            pdfDoc = pdfDoc.load(java.io.FileInputStream(pdfName));
            
            pdfStripper = PDFTextStripper();
            text = splitlines(char(pdfStripper.getText(pdfDoc)));
            
            obj = fillProps(obj, text);
           
            page = pdfDoc.getPage(0);
            pdfResource = page.getResources();
            namesArray = pdfResource.getXObjectNames().toArray();
            object = pdfResource.getXObject(namesArray(2));
            bufferedImg = object.getImage();

            H=bufferedImg.getHeight;
            W=bufferedImg.getWidth;
            
            img = uint8(zeros([H,W,3]));
            pixelsData = uint8(bufferedImg.getData.getPixels(0,0,W,H,[]));

            for i = 1 : H
                base = (i-1)*W*3+1;
                img(i,1:W,:) = deal(reshape(pixelsData(base:(base+3*W-1)),3,W)');
            end
            
            obj.image = img;
            
            pdfDoc.close();

        end
        
        function obj = clearImg(obj)

             for i = 1 : obj.initImgSize   
                for j = 1 : obj.initImgSize
                    
                    obj = changeColor(obj, 255, 255, 255, i, j);
                    obj = changeColor(obj, 0, 0, 0, i, j);
                    obj = changeColor(obj, 160, 160, 160, i, j);

                end
            end      
        end
         
        function obj = calcR(obj)
            
           myR = @(K) (1.3771 - 1) / K;
           
           obj.R1 = myR(obj.K1);
           obj.R2 = myR(obj.K2);
           obj.R1Post = myR(obj.K1Post);
           obj.R2Post = myR(obj.K2Post);
          
        end
        
        function obj = calcKPost(obj)

           obj.K1Post = obj.K1 +  distribution(obj, obj.K1Axis);
           obj.K2Post = obj.K2 +  distribution(obj, obj.K2Axis);

        end
            
        function obj = changeColor(obj, red, green, blue, i, j )

            if isequal(red, green, blue, 255)
                
                if isequal(obj.image(i, j, 1), red) && isequal(obj.image(i,j,2),green) && isequal(obj.image(i,j,3),blue)
                    obj.image(i, j, :) = obj.pixelTab(25,1,:);
                end
                
            elseif isequal(obj.image(i, j, 1), red) && isequal(obj.image(i,j,2),green) && isequal(obj.image(i,j,3),blue)
                
                y = i;
                x = j;
                
                while isequal(obj.image(y, j, 1), red) && isequal(obj.image(y,j,2),green) && isequal(obj.image(y,j,3),blue)
                    y = y + 1;
                end
                
                while isequal(obj.image(i, x, 1), red) && isequal(obj.image(i,x,2),green) && isequal(obj.image(i,x,3),blue)
                    x= x + 1;
                end
                
                tempLeft = obj.image(i, j-1, :);
                tempRight = obj.image(i, x, :);
                tempUp = obj.image(i-1, j, :);
                tempDown = obj.image(y, j, :);
               
                
                if isequal(tempLeft, tempRight, tempUp, tempDown)
                    obj.image(i,j,:) = tempLeft; 
                     
                elseif isequal(tempLeft, tempUp, tempRight)
                    obj.image(i,j,:) = tempLeft; 
                    
                elseif isequal(tempUp, tempRight, tempDown)
                    obj.image(i,j,:) = tempUp; 
                    
                elseif isequal(tempRight, tempDown, tempLeft)
                    obj.image(i,j,:) = tempLeft; 
                    
                elseif isequal(tempDown, tempLeft, tempUp)
                    obj.image(i,j,:) = tempLeft; 
                    
                elseif  isequal(tempUp,tempDown) 
                     obj.image(i, j, :) = tempUp; 
                      
                elseif  isequal(tempLeft,tempRight)
                     obj.image(i, j, :) = tempLeft; 
                     
                else
                    for x = i - 1 : i + 1
                        for y = j - 1 : j + 1
                            for k = 1 : obj.numOfColors
                                if isequal(obj.image(x, y, 1),obj.pixelTab(k, 1, 1)) && isequal(obj.image(x, y, 2),obj.pixelTab(k, 1, 2)) && isequal(obj.image(x, y, 3),obj.pixelTab(k, 1, 3)) 
                                     obj.grid(1,k) = obj.grid(1,k) + 1;
                                end
                            end 
                        end
                    end
                    
                    [val, idx] = max(obj.grid);
                    obj.grid(idx) = 0;
                    obj.image(i, j, :) = obj.pixelTab(idx,1,:);
                    obj.grid = zeros(1, obj.numOfColors);
                    
                end
            end
            
        end

        function obj = scaling(obj) 
            
            k = 0;
            
            for i = 1 : obj.initImgSize
                if isequal(obj.image(i, obj.initImgSize/2 + 1 , 1), 0) && isequal(obj.image(i,obj.initImgSize/2 + 1,2),0) && isequal(obj.image(i,obj.initImgSize/2 + 1,3),0)
                    k = i - k;
                end
            end
            
            obj.scale = obj.pupilDiameter / k;
            
        end

        function obj = processImage(obj) 
            
            obj.z = zeros(obj.initImgSize, obj.initImgSize);
            
            for i = 1 : length(obj.image)
                for j = 1 : length(obj.image)                
                    for k = 1 : obj.numOfColors
                        if isequal(obj.image(j,i,:), obj.pixelTab(k,1,:))
                            obj.z(i,j) = obj.valueTab(1,k);
                            break
                        end
                    end
                end
            end
            
        end
        
        function get3DPlot(obj)
            
            tempBorder = obj.scale*360/2;
            tempX = linspace(-tempBorder,tempBorder,360);
            tempY = linspace(-tempBorder,tempBorder,360);
            
            figure
            subplot(1,2,1);
            
            rotate(mesh(obj.x,obj.y,obj.abZ),[0 0 1], obj.K1Axis);
            xlabel('x[mm]','FontSize', 20);
            ylabel('y[mm]','FontSize', 20);
            zlabel('ablation depth[um]','FontSize', 20);
            subplot(1,2,2);
            rotate(mesh(tempX,tempY,obj.z),[0 0 1], 90);
            xlabel('x[mm]','FontSize', 20);
            ylabel('y[mm]','FontSize', 20);
            zlabel('ablation depth[um]','FontSize', 20);
            
        end
        
    end 
    
end
