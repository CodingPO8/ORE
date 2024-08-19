#include <qle/pricingengines/analytichwswaptionengine.hpp>
#include <boost/bind/bind.hpp>

#include <ostream>

namespace QuantExt {

AnalyticHwSwaptionEngine::AnalyticHwSwaptionEngine(const Array& t,
                                                    const Swaption& swaption,
                                                    const ext::shared_ptr<HwModel>&model,
                                                    const Handle<YieldTermStructure>& discountCurve)
    : GenericEngine<Swaption::arguments, Swaption::results>(), PiecewiseConstantHelper1(t), 
    model_(model), p_(model->parametrization()),
    c_(discountCurve.empty() ? p_->termStructure() : discountCurve),
    swaption(swaption){
        const Leg& origFixedLeg = swaption.underlyingSwap()->fixedLeg();
        for (const auto& cf : origFixedLeg) {
            ext::shared_ptr<FixedRateCoupon> fixedCoupon = ext::dynamic_pointer_cast<FixedRateCoupon>(cf);
            if (fixedCoupon) {
                fixedLeg_.push_back(fixedCoupon);
            }
        }
        const Leg& origFloatingLeg = swaption.underlyingSwap()->floatingLeg();
        for (const auto& cf : origFloatingLeg) {
            ext::shared_ptr<FloatingRateCoupon> coupon = ext::dynamic_pointer_cast<FloatingRateCoupon>(cf);
            if (coupon) {
                floatingLeg_.push_back(coupon); 
            }
        }
        registerWith(model);
        registerWith(c_);
}


Real AnalyticHwSwaptionEngine::s(const Time t, const Array& x) const {
    std::cout << "This is function s();" << std::endl;
    Real floatingLegNPV = 0.0;
    Real fixedLegNPV = 0.0;

    // Calculate NPV of floating leg
    for (auto const& coupon : floatingLeg_) {
        Time paymentTime = c_->timeFromReference(coupon->date());//coupon->date()就是paymentDate
        Real discountFactor = model_->discountBond(t, paymentTime, x, c_);
        floatingLegNPV += coupon->amount() * discountFactor;

        // std::cout << "Calling discountBond with t: " << t << ", T: " << paymentTime << std::endl;
        // std::cout << "discount factor: "<<discountFactor << std::endl;
        // std::cout << coupon->amount() << std::endl;
        // std::cout << std::endl;
    }
    // Calculate NPV of fixed leg
    for (auto const& coupon : fixedLeg_) {
        Time paymentTime = c_->timeFromReference(coupon->date());
        Real discountFactor = model_->discountBond(t, paymentTime, x, c_);
        fixedLegNPV += coupon->amount() * discountFactor;
        
    }
    return floatingLegNPV - fixedLegNPV;
    //return swaption.underlyingSwap()->fixedLegNPV() - swaption.underlyingSwap()->floatingLegNPV();
    
}

Real AnalyticHwSwaptionEngine::a(const Time t, const Array& x) const{
    std::cout << "This is function a();" << std::endl;
    Real annuity = 0.0;
    for (auto const& coupon : fixedLeg_){
        Time paymentTime = c_->timeFromReference(coupon->date());
        Real discountFactor = model_->discountBond(t, paymentTime, x, c_);
        annuity += coupon->accrualPeriod() * discountFactor;
    }
    return annuity; 
}

Array AnalyticHwSwaptionEngine::q(const Time t, const Array& x) const {
    std::cout << "This is function q();" << std::endl;
    Size n = p_->n();
    Array q(n, 0.0);

    Date startDate = swaption.underlyingSwap()->fixedSchedule().startDate();
    Time T0 = c_->timeFromReference(startDate);
    Date endDate = swaption.underlyingSwap()->fixedSchedule().endDate();
    Time TN = c_->timeFromReference(endDate);

    Real P_T0 = model_->discountBond(t, T0, x, c_);
    Real P_TN = model_->discountBond(t, TN, x, c_);
    
    // calculate G_j(t, T0) and G_j(t, T_N)
    for (Size j = 0; j < n; ++j) {
        Real G_j_T0 = p_->g(t, T0)[j];
        Real G_j_TN = p_->g(t, TN)[j];
        q[j] = -(P_T0 * G_j_T0 - P_TN * G_j_TN) / a(t, x);
    }

    for (Size j = 0; j < n; ++j) {
        Real sum = 0.0;
        for (const auto& coupon : fixedLeg_) {
            Time T_i = c_->timeFromReference(coupon->date());
            Real P_Ti = model_->discountBond(t, T_i, x, c_);
            Real G_j_Ti = model_->parametrization()->g(t, T_i)[j];
            sum += coupon->accrualPeriod() * P_Ti * G_j_Ti;
        }
        q[j] -= (s(t, x) * sum) / a(t, x);
    }

    return q;
}


Real AnalyticHwSwaptionEngine::v() const{
    std::cout << "This is function v();" << std::endl;
    Date startDate = swaption.underlyingSwap()->fixedSchedule().startDate();
    Time T0 = c_->timeFromReference(startDate);
    std::cout << "T0 in v(): " << T0 << std::endl;

    //PiecewiseConstantHelper1 piecewiseHelper(t_);
    
    auto integrand = [this](Real t) -> Real {
        Array x(p_->n(), 0.0); 
        QL_REQUIRE(p_->sigma_x(t).columns() == q(t,x).size(), "sigma_x and q arrays must be of the same size to perform element-wise multiplication");
        Array product = p_->sigma_x(t) * q(t,x);
        return std::pow(Norm2(product), 2);
    };

    std::cout << "print t_ in PiecewiseConstantHelper1: "<< std::endl;
    for (Time t : t_){
        std::cout << t << " ";
    }
    std::cout << std::endl;

    for (Size i = 0; i < y_->size(); ++i) {
        y_->setParam(i, inverse(integrand(t_[i])));
    }
    PiecewiseConstantHelper1::update();
    Real integral = this->int_y_sqr(T0);
    std::cout << "integral in v() is: " << integral << std::endl;

    return integral;

}

void AnalyticHwSwaptionEngine::calculate() const {
    QL_REQUIRE(arguments_.settlementType == Settlement::Physical, "cash-settled swaptions are not supported ...");
    
    Date reference = p_->termStructure()->referenceDate();
    std::cout << "reference date is: " << reference << std::endl;
    Date expiry = arguments_.exercise->dates().back();
    std::cout << "expiry date is: " << expiry << std::endl;

    if (expiry <= reference) {
        // swaption is expired, possibly generated swap is not
        // valued by this engine, so we set the npv to zero
        results_.value = 0.0;
        return;
    }
    //TODO: fix where should x come from?
    Array x(p_->n(), 0.0);

    NormalDistribution pdf;
    CumulativeNormalDistribution cdf;

    const auto& underlyingSwap = swaption.underlyingSwap();
    Real strikePrice = underlyingSwap->fixedRate();

    Real s0 = s(0, x);
    std::cout << "s0 is: " << s0 << std::endl;
    Real a0 = a(0, x);
    std::cout << "a0 is: " << a0 << std::endl;
    Real c = strikePrice; 
    Real d = (s0 - c) / sqrt(v());
    std::cout << "d is: " << std::endl;
    results_.value = a0 * ((s0 - c) * cdf(d) + sqrt(v()) * pdf(d));
    std::cout << "This is the final result!: " << results_.value << std::endl; 
}

}