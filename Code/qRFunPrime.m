function qqRPrime = qRFunPrime(gamma,r)
    % Analytical formula for calculating q_r'(gamma,r) - Obtained by symbolic solver
    % Require gamma < 1, which should always occur when a dielectric exists. This is because gamma = w0 / (d0 + di/epsr), and the maximum value of w0 is d0.
    % 0<= r <= 1, qr(gamma,r=1) = q(gamma)
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
        r (1,1) {mustBeNumeric,mustBeLessThanOrEqual(r,1),mustBeGreaterThanOrEqual(r,0)} = 1 %when r = 1, qR(gamma,r) = q(gamma)
    end
    
    %analytical derivative of qR(gamma,r)
    term1 = 2* r *sqrt(gamma) .* (1 + gamma -r^2 * gamma) ./ ( (gamma-1) .* (-1 + (r^2-1)^2 *gamma ));
    term2 = (2*sqrt(gamma) - 1) .* sqrt(sqrt(gamma)-gamma) .* atan(r*sqrt(gamma) ./ sqrt(sqrt(gamma)-gamma)) ./ (sqrt(gamma)-1).^2;
    term3 = -(1+2*sqrt(gamma)) .* sqrt(sqrt(gamma)+gamma) .* atanh(r*sqrt(gamma) ./(sqrt(sqrt(gamma)+gamma))) ./ (1+sqrt(gamma)).^2;

    qqRPrime = 1./(4* gamma.^ (3/2) ) .* (term1+term2 + term3);

    % The critical point to evaluate is @ gamma = 0.
    try
        ZeroIndex = find(abs(gamma) < eps); %to prevent floating point errors where gamma is extremely close to zero.
        % eps is 2.2 e-16 so no inputs should ever be below this outside of 0.
    catch
        disp('problem with eps when evaluating qFun')
        ZeroIndex = find (gamma == 0); %in case there is a problem with eps, which should not happen since functions have their own workspace.
    end
    
    if ~isempty(ZeroIndex)
        qqRPrime(ZeroIndex) = 2 *r - 4*r^3/3 + 2 * r^5 /5;
        %note: qqPrime(0) = 2 *r - 4*r^3/3 + 2 * r^5 /5 if you evaluate the limit, but MATLAB just evaluates 0/0 = NaN.
    end
end %end function