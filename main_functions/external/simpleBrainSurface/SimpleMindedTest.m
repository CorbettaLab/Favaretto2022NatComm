classdef SimpleMindedTest < matlab.unittest.TestCase
    
    properties
       TestFigure 
    end
    
    methods(TestMethodSetup)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % newFigure
        % Setting up a new figure for testing
        % and assinging to testCase property
        % INPUT:
        % testCase
        % OUTPUT:
        % None.
        % SIDEEFFECTS:
        % new figure is created and handle stored in testCase.TestFigure
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function newFigure(testCase)
           testCase.TestFigure = figure; 
        end
    end
    
    methods(TestMethodTeardown)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %closeFigure
        % Deletes the set up test figure.
        % INPUT:
        % OUTPUT:
        % SIDEEFFECTS:
        % figure is deleted.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function closeFigure(testCase)
            close(testCase.TestFigure);
        end
    end
   
    methods(Test)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % checkBrain
        % Checks simply if the function draws the right
        % number of triangles based on the number of vertices
        % and patches.
        % INPUT:
        % testCase
        % OUTPUT:
        % None.
        % SIDEEFFECTS:
        % brain surface figure is generated
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function checkBrain(testCase)
           h = simpleBrainSurface();
           testCase.verifyEqual(size(get(h,'Vertices')),[81924,3]);
           testCase.verifyEqual(size(get(h,'Faces')),[163840,3]);
        end
        
        function checkDefaultShading(testCase)
           h = simpleBrainSurface();
           vertexcdata = get(h,'facevertexcdata');
           testCase.verifyTrue(all(all(vertexcdata(:,1)==vertexcdata(:,2) & ...
               vertexcdata(:,1)==vertexcdata(:,3))));
           testCase.check_min_max(0.1,0.7,vertexcdata(:,1));
        end
        
        function checkSpecsObject(testCase)
           specs_object.range = [0.05 0.8]; 
           h = simpleBrainSurface(specs_object);
           vertexcdata = get(h,'facevertexcdata');
           testCase.verifyTrue(all(all(vertexcdata(:,1)==vertexcdata(:,2) & ...
               vertexcdata(:,1)==vertexcdata(:,3))));
           testCase.check_min_max(0.05,0.8,vertexcdata(:,1));
        end
        
    end
    
    methods(Access = private)
        
        function check_min_max(testCase, min_val, max_val, data)
           testCase.verifyEqual(min(data),min_val,'AbsTol',0.001);
           testCase.verifyEqual(max(data),max_val,'AbsTol',0.001);
        end
        
    end % private methods
    
end