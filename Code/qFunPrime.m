function qqPrime = qFunPrime(gamma)
    % Formula for calculating q'(gamma) - eqn 15 in manuscript.
    % Require gamma < 1, which should always occur when a dielectric exists. TThis is because gamma = w0 / (d0 + di/epsr), and the maximum value of w0 is d0.
    %
    % Copyright 2024 Eric Dew
    %
    % Licensed under the Apache License, Version 2.0 (the "License");
    % you may not use this file except in compliance with the License.
    % You may obtain a copy of the License at
    %
    %     https://www.apache.org/licenses/LICENSE-2.0
    %
    % Unless required by applicable law or agreed to in writing, software
    % distributed under the License is distributed on an "AS IS" BASIS,
    % WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    % See the License for the specific language governing permissions and
    % limitations under the License.
    arguments
        gamma (1,:) {mustBeNumeric,mustBeLessThan(gamma,1)}       
    end

    %analytical derivative of q(gamma)
    qqPrime =  (1 ./ (4 * (gamma - 1).^2 .* gamma .^ (3/2))) .* ...
        ( sqrt( sqrt(gamma) -gamma) .* (2 * gamma .^ (3/2) + 3*gamma -1) .* ...
        atan(  sqrt(gamma)./ (sqrt(sqrt(gamma) -gamma)) ) -...
        sqrt(gamma + sqrt(gamma)) .* (2* gamma.^(3/2) - 3 * gamma + 1) .* ...
        atanh( sqrt(gamma)./ (sqrt(sqrt(gamma) + gamma)) ) -...
        2*sqrt(gamma) .* (gamma-1)); %analytical derivative    

    % The critical point to evaluate is @ gamma = 0.
    try
        ZeroIndex = find(abs(gamma) < eps); %to prevent floating point errors where gamma is extremely close to zero.
        % eps is 2.2 e-16 so no inputs should ever be below this outside of 0.
    catch
        disp('problem with eps when evaluating qFun')
        ZeroIndex = find (gamma == 0); %in case there is a problem with eps, which should not happen since functions have their own workspace.
    end
    
    if ~isempty(ZeroIndex)
        qqPrime(ZeroIndex) = 16/15;
        %note: qqPrime(0) = 16/15 if you evaluate the limit, but MATLAB just evaluates 0/0 = NaN.
    end
end %end function