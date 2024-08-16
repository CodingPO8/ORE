#include <qle/pricingengines/analyticnfhwswaptionengine.hpp>
#include <qle/models/hwmodel.hpp>
#include <qle/models/hwparametrization.hpp>
#include <ql/math/integrals/integral.hpp>
#include <ql/math/matrix.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/exercise.hpp>

#include <boost/bind.hpp>

#include <cmath>

namespace QuantExt {

    AnalyticNfhwSwaptionEngine::AnalyticNfhwSwaptionEngine(const boost::shared_ptr<HwModel>& model,
        const Handle<YieldTermStructure>& discountCurve)
        : GenericEngine<Swaption::arguments, Swaption::results>(), p_(model->parametrization()),
        c_(discountCurve.empty() ? p_->termStructure() : discountCurve){
        registerWith(model);
        registerWith(c_);
    }

    AnalyticNfhwSwaptionEngine::AnalyticNfhwSwaptionEngine(const boost::shared_ptr<IrHwParametrization>& parametrization,
        const Handle<YieldTermStructure>& discountCurve)
        : GenericEngine<Swaption::arguments, Swaption::results>(), p_(parametrization),
        c_(discountCurve.empty() ? p_->termStructure() : discountCurve){
        registerWith(c_);
    }

    // myFunction will be used to calculate the sum of squares part in the formula
    Real AnalyticNfhwSwaptionEngine::myFunction(const Time t, const Real x, const Handle<YieldTermStructure>& c_) const{
        // create a hwmodel based on the p_ for the use of discountBond
        Date T_0 = floatingLeg_[0]->date();
        Date T_N = floatingLeg_[floatingLeg_.size()-1]->date();
        Real P_0 = discountBond(t, T_0, x, c_);
        Real P_N = discountBond(t, T_N, x, c_);
        Real A_t = 0.0;
        for (Size i = 0; i < floatingLeg_.size()-1; ++i) {
            Real discountFactor = discountBond(t, floatingLeg_[i+1]->date(), x, c_);
            Real Tau = floatingLeg_[i]->accrualPeriod();
            A_t += Tau * discountFactor;
        }
        Real S_t = (P_0 - P_N) / A_t;
        Array SumTemp; 
        for (Size i = 0; i < floatingLeg_.size()-1; ++i) {
            Real discountFactor =discountBond(t, floatingLeg_[i+1]->date(), x, c_);
            Array gt_temp = p_->g(t, floatingLeg_[i+1]->date()); 
            Real Tau = floatingLeg_[i]->accrualPeriod();
            SumTemp = SumTemp + Tau * gt_temp * discountFactor;
        }
        Array gt_0 = p_->g(t, floatingLeg_[0]->date());
        Array gt_N = p_->g(t, floatingLeg_[floatingLeg_.size()-1]->date());
        // Make a matrix from the array 
        // Matrix gt_0_as_matrix(1,gt_0.size(),gt_0.begin(),gt_0.end())
        // make sure the sigma_x(t) is also a matrix
        Matrix product_result = transpose(-(P_0 * gt_0 - P_N * gt_N) / A_t - S_t * SumTemp / A_t) * transpose(p_->sigma_x(t));
        // get product_result(0,0) 
        double sumOfSquares = 0.0;
        for (double value : product_result) {
            sumOfSquares += value * value;
        }
        return std::sqrt(sumOfSquares);
    }

    void AnalyticNfhwSwaptionEngine::calculate() const {
        //parametrization

        QL_REQUIRE(arguments_.settlementType == Settlement::Physical, "cash-settled swaptions are not supported ...");

        Date reference = p_->termStructure()->referenceDate();// give us the date 
        //a function under termStructure(), timefromreference() converts a date to a time

        Date expiry = arguments_.exercise->dates().back();
        Time expiry2 = p_->termStructure()->timeFromReference(expiry);

        if (expiry <= reference) {
            // swaption is expired, possibly generated swap is not
            // valued by this engine, so we set the npv to zero
            results_.value = 0.0;
            return;
        }
        // 1. Calculate A(0) and S(0)
        Real x = 0.0; //approximation
        Real A_0 = 0.0;
        for (Size i = 0; i < floatingLeg_.size()-1; ++i) {
            Real discountFactor =discountBond(0.0, floatingLeg_[i+1]->date(), x, c_);
            Real Tau = floatingLeg_[i]->accrualPeriod();
            A_0 += Tau * discountFactor;
        }

        Real S_0 = (discountBond(0.0, floatingLeg_[0]->date(), x, c_) - discountBond(0.0, floatingLeg_[floatingLeg_.size() - 1]->date(), x, c_)) / A_0;

        // 2. Calculate the swaption price
        // upper limit should be the "expiry", change this to a real number
        // for the integrator function, we can use the segment integral (segmentintegral.hpp)
        
        Real v = integrator(myFunction, 0.0, expiry2);
        Real d = (S_0 - fixedLeg_[0]->rate()) / sqrt(v);
        // Create an object of normal distribution to use cdf, pdf function
        NormalDistribution pdf;
        CumulativeNormalDistribution cdf;
        Real cdf_d = cdf(d);
        Real pdf_d = pdf(d);
        Real sum = (A_0 * ((S_0 - fixedLeg_[0]->rate()) * cdf_d + sqrt(v) * pdf_d));
        results_.value = sum;
        return; // swaption price
    }

} // namespace QuantExt
