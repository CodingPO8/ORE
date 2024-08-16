#pragma once

#ifndef quantext_analytic_nfhw_swaptionengine_hpp
#define quantext_analytic_nfhw_swaptionengine_hpp

#include <ql/instruments/swaption.hpp>
#include <ql/cashflows/floatingratecoupon.hpp>
#include <ql/cashflows/fixedratecoupon.hpp>
#include <qle/models/hwmodel.hpp>
#include <qle/models/hwparametrization.hpp>

namespace QuantExt {

    class AnalyticNfhwSwaptionEngine : public GenericEngine<Swaption::arguments, Swaption::results> {

    public:
        //enum class StateApproximation { Zero, FirstOrder };

        // Constructor accepting HwModel
        AnalyticNfhwSwaptionEngine(const boost::shared_ptr<HwModel>& model,
            const Handle<YieldTermStructure>& discountCurve = Handle<YieldTermStructure>());


        // Constructor accepting hwparametrization
        AnalyticNfhwSwaptionEngine(const boost::shared_ptr<IrHwParametrization>& parametrization,
            const Handle<YieldTermStructure>& discountCurve = Handle<YieldTermStructure>());


        void calculate() const override;

    private:
        // Private member variables
        const boost::shared_ptr<IrHwParametrization> p_;
        const Handle<YieldTermStructure> c_;
        mutable std::vector<boost::shared_ptr<FixedRateCoupon>> fixedLeg_;
        mutable std::vector<boost::shared_ptr<FloatingRateCoupon>> floatingLeg_;
        Real myFunction(const Real t, const Real x, const Handle<YieldTermStructure>& c_) const;
    };

} // namespace QuantExt

#endif

