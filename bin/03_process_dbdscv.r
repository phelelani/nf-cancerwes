<!DOCTYPE html><html lang="en-US" class="" data-primer><head><link href="https://a.slack-edge.com/d5fba4c/marketing/style/onetrust/onetrust_banner.css" rel="stylesheet" type="text/css" onload="window._cdn ? _cdn.ok(this, arguments) : null" onerror="window._cdn ? _cdn.failed(this, arguments) : null" crossorigin="anonymous"><link href="https://a.slack-edge.com/e06451a/style/libs/lato-2-compressed.css" rel="stylesheet" type="text/css" onload="window._cdn ? _cdn.ok(this, arguments) : null" onerror="window._cdn ? _cdn.failed(this, arguments) : null" crossorigin="anonymous"><link href="https://a.slack-edge.com/css/v5/style/_generic.typography.larsseit.85ad0e0bbe61bdbf62bdd9efa15a921e01033c37.css" rel="stylesheet" type="text/css" onload="window._cdn ? _cdn.ok(this, arguments) : null" onerror="window._cdn ? _cdn.failed(this, arguments) : null" crossorigin="anonymous"><link rel="canonical" href="https://slack.com">

<link rel="alternate" hreflang="en-us" href="https://slack.com">

<link rel="alternate" hreflang="en-au" href="https://slack.com/intl/en-au">

<link rel="alternate" hreflang="en-gb" href="https://slack.com/intl/en-gb">

<link rel="alternate" hreflang="en-in" href="https://slack.com/intl/en-in">

<link rel="alternate" hreflang="fr-ca" href="https://slack.com/intl/fr-ca">

<link rel="alternate" hreflang="fr-fr" href="https://slack.com/intl/fr-fr">

<link rel="alternate" hreflang="de-de" href="https://slack.com/intl/de-de">

<link rel="alternate" hreflang="es-es" href="https://slack.com/intl/es-es">

<link rel="alternate" hreflang="es" href="https://slack.com/intl/es-la">

<link rel="alternate" hreflang="ja-jp" href="https://slack.com/intl/ja-jp">

<link rel="alternate" hreflang="pt-br" href="https://slack.com/intl/pt-br">

<link rel="alternate" hreflang="ko-kr" href="https://slack.com/intl/ko-kr">

<link rel="alternate" hreflang="it-it" href="https://slack.com/intl/it-it">

<link rel="alternate" hreflang="zh-cn" href="https://slack.com/intl/zh-cn">

<link rel="alternate" hreflang="zh-tw" href="https://slack.com/intl/zh-tw">

<link rel="alternate" hreflang="x-default" href="https://slack.com">

<script>window.ts_endpoint_url = "https:\/\/slack.com\/beacon\/timing";(function(e) {
	var n=Date.now?Date.now():+new Date,r=e.performance||{},t=[],a={},i=function(e,n){for(var r=0,a=t.length,i=[];a>r;r++)t[r][e]==n&&i.push(t[r]);return i},o=function(e,n){for(var r,a=t.length;a--;)r=t[a],r.entryType!=e||void 0!==n&&r.name!=n||t.splice(a,1)};r.now||(r.now=r.webkitNow||r.mozNow||r.msNow||function(){return(Date.now?Date.now():+new Date)-n}),r.mark||(r.mark=r.webkitMark||function(e){var n={name:e,entryType:"mark",startTime:r.now(),duration:0};t.push(n),a[e]=n}),r.measure||(r.measure=r.webkitMeasure||function(e,n,r){n=a[n].startTime,r=a[r].startTime,t.push({name:e,entryType:"measure",startTime:n,duration:r-n})}),r.getEntriesByType||(r.getEntriesByType=r.webkitGetEntriesByType||function(e){return i("entryType",e)}),r.getEntriesByName||(r.getEntriesByName=r.webkitGetEntriesByName||function(e){return i("name",e)}),r.clearMarks||(r.clearMarks=r.webkitClearMarks||function(e){o("mark",e)}),r.clearMeasures||(r.clearMeasures=r.webkitClearMeasures||function(e){o("measure",e)}),e.performance=r,"function"==typeof define&&(define.amd||define.ajs)&&define("performance",[],function(){return r}) // eslint-disable-line
})(window);</script><script>

(function () {
	
	window.TSMark = function (mark_label) {
		if (!window.performance || !window.performance.mark) return;
		performance.mark(mark_label);
	};
	window.TSMark('start_load');

	
	window.TSMeasureAndBeacon = function (measure_label, start_mark_label) {
		if (!window.performance || !window.performance.mark || !window.performance.measure) {
			return;
		}

		performance.mark(start_mark_label + '_end');

		try {
			performance.measure(measure_label, start_mark_label, start_mark_label + '_end');
			window.TSBeacon(measure_label, performance.getEntriesByName(measure_label)[0].duration);
		} catch (e) {
			
		}
	};

	
	if ('sendBeacon' in navigator) {
		window.TSBeacon = function (label, value) {
			var endpoint_url = window.ts_endpoint_url || 'https://slack.com/beacon/timing';
			navigator.sendBeacon(
				endpoint_url + '?data=' + encodeURIComponent(label + ':' + value),
				'',
			);
		};
	} else {
		window.TSBeacon = function (label, value) {
			var endpoint_url = window.ts_endpoint_url || 'https://slack.com/beacon/timing';
			new Image().src = endpoint_url + '?data=' + encodeURIComponent(label + ':' + value);
		};
	}
})();
</script><script>window.TSMark('step_load');</script><script>
(function () {
	function throttle(callback, intervalMs) {
		var wait = false;

		return function () {
			if (!wait) {
				callback.apply(null, arguments);
				wait = true;
				setTimeout(function () {
					wait = false;
				}, intervalMs);
			}
		};
	}

	function getGenericLogger() {
		return {
			warn: (msg) => {
				
				if (window.console && console.warn) return console.warn(msg);
			},
			error: (msg) => {
				if (!msg) return;

				if (window.TSBeacon) return window.TSBeacon(msg, 1);

				
				if (window.console && console.warn) return console.warn(msg);
			},
		};
	}

	function globalErrorHandler(evt) {
		if (!evt) return;

		
		var details = '';

		var node = evt.srcElement || evt.target;

		var genericLogger = getGenericLogger();

		
		
		
		
		if (node && node.nodeName) {
			var nodeSrc = '';
			if (node.nodeName === 'SCRIPT') {
				nodeSrc = node.src || 'unknown';
				var errorType = evt.type || 'error';

				
				details = `[${errorType}] from script at ${nodeSrc} (failed to load?)`;
			} else if (node.nodeName === 'IMG') {
				nodeSrc = node.src || node.currentSrc;

				genericLogger.warn(`<img> fired error with url = ${nodeSrc}`);
				return;
			}
		}

		
		if (evt.error && evt.error.stack) {
			details += ` ${evt.error.stack}`;
		}

		if (evt.filename) {
			
			var fileName = evt.filename;
			var lineNo = evt.lineno;
			var colNo = evt.colno;

			details += ` from ${fileName}`;

			
			if (lineNo) {
				details += ` @ line ${lineNo}, col ${colNo}`;
			}
		}

		var msg;

		
		if (evt.error && evt.error.stack) {
			
			msg = details;
		} else {
			
			msg = `${evt.message || ''} ${details}`;
		}

		
		if (msg && msg.replace) msg = msg.replace(/\s+/g, ' ').trim();

		if (!msg || !msg.length) {
			if (node) {
				var nodeName = node.nodeName || 'unknown';

				
				
				
				if (nodeName === 'VIDEO') {
					return;
				}

				msg = `error event from node of ${nodeName}, no message provided?`;
			} else {
				msg = 'error event fired, no relevant message or node found';
			}
		}

		var logPrefix = 'ERROR caught in js/inline/register_global_error_handler';

		msg = `${logPrefix} - ${msg}`;

		genericLogger.error(msg);
	}

	
	
	
	var capture = true;

	
	var throttledHandler = throttle(globalErrorHandler, 10000);

	window.addEventListener('error', throttledHandler, capture);
})();
</script><script type="text/javascript" crossorigin="anonymous" src="https://a.slack-edge.com/bv1-13/manifest.5104750ec1047909eb90.primer.min.js" onload="window._cdn ? _cdn.ok(this, arguments) : null" onerror="window._cdn ? _cdn.failed(this, arguments) : null"></script><noscript><meta http-equiv="refresh" content="0; URL=/?redir=%2Ffiles%2FUBTGAQJP2%2FF0755HYPHGE%2F03_process_dbdscv.r%3Fu%3DUBTGAQJP2%26file_id%3DF0755HYPHGE%26name%3D03_process_dbdscv.r&amp;nojsmode=1"></noscript><script type="text/javascript">var safe_hosts = ['app.optimizely.com', 'tinyspeck.dev.slack.com'];

if (self !== top && safe_hosts.indexOf(top.location.host) === -1) {
	window.document.write(
		'\u003Cstyle>body * {display:none !important;}\u003C/style>\u003Ca href="#" onclick=' +
			'"top.location.href=window.location.href" style="display:block !important;padding:10px">Go to Slack.com\u003C/a>'
	);
}

(function() {
	var timer;
	if (self !== top && safe_hosts.indexOf(top.location.host) === -1) {
		timer = window.setInterval(function() {
			if (window) {
				try {
					var pageEl = document.getElementById('page');
					var clientEl = document.getElementById('client-ui');
					var sectionEls = document.querySelectorAll('nav, header, section');

					pageEl.parentNode.removeChild(pageEl);
					clientEl.parentNode.removeChild(clientEl);
					for (var i = 0; i < sectionEls.length; i++) {
						sectionEls[i].parentNode.removeChild(sectionEls[i]);
					}
					window.TS = null;
					window.TD = null;
					window.clearInterval(timer);
				} catch (e) {}
			}
		}, 200);
	}
})();</script><meta name="facebook-domain-verification" content="chiwsajpoybn2cnqyj9w8mvrey56m0"><script type="text/javascript">
document.addEventListener("DOMContentLoaded", function(e) {
	var gtmDataLayer = window.dataLayer || [];
	var gtmTags = document.querySelectorAll('*[data-gtm-click]');
	var gtmClickHandler = function(c) {
		var gtm_events = this.getAttribute('data-gtm-click');
		if (!gtm_events) return;
		var gtm_events_arr = gtm_events.split(",");
		for(var e=0; e < gtm_events_arr.length; e++) {
			var ev = gtm_events_arr[e].trim();
			gtmDataLayer.push({ 'event': ev });
		}
	};
	for(var g=0; g < gtmTags.length; g++){
		var elem = gtmTags[g];
		elem.addEventListener('click', gtmClickHandler);
	}
});
</script><script type="text/javascript">
(function(e,c,b,f,d,g,a){e.SlackBeaconObject=d;
e[d]=e[d]||function(){(e[d].q=e[d].q||[]).push([1*new Date(),arguments])};
e[d].l=1*new Date();g=c.createElement(b);a=c.getElementsByTagName(b)[0];
g.async=1;g.src=f;a.parentNode.insertBefore(g,a)
})(window,document,"script","https://a.slack-edge.com/bv1-13/slack_beacon.b9cabcadb13e05735d6c.min.js","sb");
window.sb('set', 'token', '3307f436963e02d4f9eb85ce5159744c');
window.sb('track', 'pageview');
</script><script src="https://cdn.cookielaw.org/scripttemplates/otSDKStub.js" data-document-language="true" data-domain-script="3bcd90cf-1e32-46d7-adbd-634f66b65b7d"></script><script>window.OneTrustLoaded = true;</script><script>
window.dataLayer = window.dataLayer || [];

function maybeShowAnnouncementBanner() {
	const bottomBannerEl = document.querySelector('.c-announcement-banner-bottom');
	if (bottomBannerEl !== null) {
		bottomBannerEl.classList.remove('c-announcement-banner-bottom-invisible');
	}
}

function OptanonWrapper() {
	gtag('consent', "update", {"ad_storage":"denied","ad_user_data":"denied","ad_personalization":"denied","personalization_storage":"denied","analytics_storage":"denied","functionality_storage":"denied","security_storage":"denied"});
	window.dataLayer.push({'event': 'OneTrustReady'});
	if (!Optanon.GetDomainData().ShowAlertNotice || false) {
		maybeShowAnnouncementBanner();
	} else {
		document.querySelector('#onetrust-accept-btn-handler').focus()
	}
	Optanon.OnConsentChanged(function() {
		maybeShowAnnouncementBanner();
	});
}</script><script type="text/javascript">var TS_last_log_date = null;
var TSMakeLogDate = function() {
	var date = new Date();

	var y = date.getFullYear();
	var mo = date.getMonth()+1;
	var d = date.getDate();

	var time = {
	  h: date.getHours(),
	  mi: date.getMinutes(),
	  s: date.getSeconds(),
	  ms: date.getMilliseconds()
	};

	Object.keys(time).map(function(moment, index) {
		if (moment == 'ms') {
			if (time[moment] < 10) {
				time[moment] = time[moment]+'00';
			} else if (time[moment] < 100) {
				time[moment] = time[moment]+'0';
			}
		} else if (time[moment] < 10) {
			time[moment] = '0' + time[moment];
		}
	});

	var str = y + '/' + mo + '/' + d + ' ' + time.h + ':' + time.mi + ':' + time.s + '.' + time.ms;
	if (TS_last_log_date) {
		var diff = date-TS_last_log_date;
		//str+= ' ('+diff+'ms)';
	}
	TS_last_log_date = date;
	return str+' ';
}

var parseDeepLinkRequest = function(code) {
	var m = code.match(/"id":"([CDG][A-Z0-9]{8,})"/);
	var id = m ? m[1] : null;

	m = code.match(/"team":"(T[A-Z0-9]{8,})"/);
	var team = m ? m[1] : null;

	m = code.match(/"message":"([0-9]+\.[0-9]+)"/);
	var message = m ? m[1] : null;

	return { id: id, team: team, message: message };
}

if ('rendererEvalAsync' in window) {
	var origRendererEvalAsync = window.rendererEvalAsync;
	window.rendererEvalAsync = function(blob) {
		try {
			var data = JSON.parse(decodeURIComponent(atob(blob)));
			if (data.code.match(/handleDeepLink/)) {
				var request = parseDeepLinkRequest(data.code);
				if (!request.id || !request.team || !request.message) return;

				request.cmd = 'channel';
				TSSSB.handleDeepLinkWithArgs(JSON.stringify(request));
				return;
			} else {
				origRendererEvalAsync(blob);
			}
		} catch (e) {
		}
	}
}</script><script type="text/javascript">var TSSSB = {
	call: function() {
		return false;
	}
};</script><title>Slack</title><meta name="referrer" content="no-referrer"><meta name="superfish" content="nofish"><meta name="author" content="Slack"><meta name="description" content=""><meta name="keywords" content=""></head><body class="full_height"><div id="page_contents"><div id="props_node" data-props="{&quot;loggedInTeams&quot;:[],&quot;entryPoint&quot;:&quot;&quot;,&quot;isPaidTeam&quot;:false,&quot;teamName&quot;:&quot;Sydney Brenner Institute for Molecular Bioscience&quot;,&quot;teamDomain&quot;:&quot;sbimb&quot;,&quot;encodedTeamId&quot;:&quot;T19L623NY&quot;,&quot;emailJustSent&quot;:false,&quot;shouldRedirect&quot;:true,&quot;redirectURL&quot;:&quot;\/files\/UBTGAQJP2\/F0755HYPHGE\/03_process_dbdscv.r?u=UBTGAQJP2&amp;file_id=F0755HYPHGE&amp;name=03_process_dbdscv.r&quot;,&quot;redirectQs&quot;:&quot;\/?redir=%2Ffiles%2FUBTGAQJP2%2FF0755HYPHGE%2F03_process_dbdscv.r%3Fu%3DUBTGAQJP2%26file_id%3DF0755HYPHGE%26name%3D03_process_dbdscv.r&quot;,&quot;remember&quot;:false,&quot;hasRemember&quot;:true,&quot;noSSO&quot;:false,&quot;crumbValue&quot;:&quot;s-1716878842-2aa65e7c3a83156162c99903c3067bc458997b4f95857f01939a596c00c2fcb7-\u2603&quot;,&quot;isSSB&quot;:false,&quot;isSSBSignIn&quot;:false,&quot;isSSBSonicSignIn&quot;:false,&quot;SSBVersion&quot;:&quot;&quot;,&quot;hasEmailError&quot;:false,&quot;emailValue&quot;:&quot;&quot;,&quot;hasPasswordError&quot;:false,&quot;isMobileBrowser&quot;:false,&quot;hasAuthReloginFlow&quot;:false,&quot;hasRateLimit&quot;:false,&quot;recaptchaSitekey&quot;:&quot;6LcQQiYUAAAAADxJHrihACqD5wf3lksm9jbnRY5k&quot;,&quot;hasSmsRateLimit&quot;:false,&quot;forgotPasswordLink&quot;:&quot;\/forgot&quot;,&quot;showSignupEmailLink&quot;:true,&quot;getStartedLink&quot;:&quot;https:\/\/slack.com\/get-started#\/find&quot;,&quot;isSSOAuthMode&quot;:false,&quot;isNormalAuthMode&quot;:true,&quot;signinUrl&quot;:&quot;https:\/\/slack.com\/signin&quot;,&quot;signinFindUrl&quot;:&quot;https:\/\/slack.com\/signin\/find&quot;,&quot;ssbRelogin&quot;:&quot;&quot;,&quot;isLoggedOutSSOMaybeRequired&quot;:false,&quot;isLoggedOutRedirect&quot;:true,&quot;teamAuthMode&quot;:null,&quot;authModeGoogle&quot;:&quot;google&quot;,&quot;samlProviderLabel&quot;:null,&quot;errorSource&quot;:&quot;&quot;,&quot;errorMissing&quot;:false,&quot;errorNoUser&quot;:false,&quot;errorDeleted&quot;:false,&quot;errorPassword&quot;:false,&quot;errorSSONoOwner&quot;:false,&quot;errorSSONonRA&quot;:false,&quot;errorTwoFactorWrong&quot;:false,&quot;errorSmsRateLimit&quot;:false,&quot;errorConfirmed&quot;:false,&quot;errorNoPassword&quot;:false,&quot;errorTwoFactorState&quot;:false,&quot;hasEmailOnDomain&quot;:false,&quot;truncatedEmailDomains&quot;:null,&quot;truncatedEmailDomainsCount&quot;:0,&quot;formattedEmailDomains&quot;:&quot;&quot;,&quot;isReloginFlow&quot;:false,&quot;downloadLinkSigninCTA&quot;:{&quot;linkUrl&quot;:&quot;\/get-started#\/create&quot;,&quot;isDownload&quot;:false},&quot;twoFactorRequired&quot;:false,&quot;usingBackup&quot;:null,&quot;twoFactorType&quot;:null,&quot;twoFactorBackupType&quot;:null,&quot;twoFactorRequiredML&quot;:null,&quot;authcodeInputType&quot;:&quot;text&quot;,&quot;backupPhoneNumber&quot;:null,&quot;forgotPasswordError&quot;:&quot;&quot;,&quot;resetLinkSent&quot;:false,&quot;userOauth&quot;:{&quot;apple&quot;:{&quot;enabled&quot;:true,&quot;redirect_url&quot;:&quot;https:\/\/sbimb.slack.com\/workspace-signin\/oauth\/apple\/start?redir=%2Ffiles%2FUBTGAQJP2%2FF0755HYPHGE%2F03_process_dbdscv.r%3Fu%3DUBTGAQJP2%26file_id%3DF0755HYPHGE%26name%3D03_process_dbdscv.r&amp;is_ssb_browser_signin=&quot;},&quot;google&quot;:{&quot;enabled&quot;:true,&quot;redirect_url&quot;:&quot;https:\/\/sbimb.slack.com\/workspace-signin\/oauth\/google\/start?redir=%2Ffiles%2FUBTGAQJP2%2FF0755HYPHGE%2F03_process_dbdscv.r%3Fu%3DUBTGAQJP2%26file_id%3DF0755HYPHGE%26name%3D03_process_dbdscv.r&amp;is_ssb_browser_signin=&quot;}},&quot;isUrgentBannerExpOn&quot;:false,&quot;isWarningBannerExpOn&quot;:true,&quot;signInWithPassword&quot;:false}"></div></div><script type="text/javascript">
/**
 * A placeholder function that the build script uses to
 * replace file paths with their CDN versions.
 *
 * @param {String} file_path - File path
 * @returns {String}
 */
function vvv(file_path) {
		 var vvv_warning = 'You cannot use vvv on dynamic values. Please make sure you only pass in static file paths.'; if (window.TS && window.TS.warn) { window.TS.warn(vvv_warning); } else { console.warn(vvv_warning); } 
	return file_path;
}

var cdn_url = "https:\/\/a.slack-edge.com";
var vvv_abs_url = "https:\/\/slack.com\/";
var inc_js_setup_data = {
	emoji_sheets: {
		apple: 'https://a.slack-edge.com/80588/img/emoji_2017_12_06/sheet_apple_64_indexed_256.png',
		google: 'https://a.slack-edge.com/80588/img/emoji_2017_12_06/sheet_google_64_indexed_256.png',
	},
};
</script><script nonce="" type="text/javascript">	// common boot_data
	var boot_data = {"cdn":{"edges":["https:\/\/a.slack-edge.com\/","https:\/\/b.slack-edge.com\/","https:\/\/a.slack-edge.com\/"],"avatars":"https:\/\/ca.slack-edge.com\/","downloads":"https:\/\/downloads.slack-edge.com\/","files":"https:\/\/slack-files.com\/"},"feature_builder_story_step":false,"feature_tinyspeck":false,"feature_olug_remove_required_workspace_setting":false,"feature_deprecate_get_member_by_name":false,"feature_file_threads":true,"feature_broadcast_indicator":true,"feature_sonic_emoji":true,"feature_attachments_inline":false,"feature_desktop_symptom_events":true,"feature_gdpr_user_join_tos":true,"feature_user_invite_tos_april_2018":true,"feature_channel_mgmt_message_count":false,"feature_channel_exports":false,"feature_allow_intra_word_formatting":true,"feature_slim_scrollbar":false,"feature_edge_upload_proxy_check":false,"feature_app_views_reminders":true,"feature_reminders_org_shard":false,"feature_blocks_reminders_list":false,"feature_message_blocks":false,"feature_set_tz_automatically":true,"feature_attachments_v2":true,"feature_beacon_js_errors":false,"feature_user_app_disable_speed_bump":true,"feature_apps_manage_permissions_scope_changes":true,"feature_ia_member_profile":true,"feature_desktop_reload_on_generic_error":true,"feature_desktop_extend_app_menu":true,"feature_desktop_restart_service_worker":false,"feature_wta_stop_creation":true,"feature_admin_email_change_confirm":false,"feature_improved_email_rendering":true,"feature_recent_desktop_files":true,"feature_cea_allowlist_changes":false,"feature_cea_channel_management":true,"feature_cea_admin_controls":false,"feature_cea_allowlist_changes_plus":false,"feature_ia_layout":true,"feature_threaded_call_block":true,"feature_enterprise_mobile_device_check":true,"feature_trace_jq_init":true,"feature_seven_days_email_update":true,"feature_indonesia_tax_change_notification":false,"feature_indonesia_tax_assessment":false,"feature_channel_sections":true,"feature_show_email_forwarded_by":false,"feature_mpdm_audience_expansion":true,"feature_remove_email_preview_link":true,"feature_desktop_enable_tslog":false,"feature_email_determine_charset":true,"feature_no_deprecation_in_updater":false,"feature_pea_domain_allowlist":true,"feature_composer_auth_admin":true,"experiment_assignments":{"cust_acq_homepage_above_the_fold_cta":{"experiment_id":"6887647122804","type":"visitor","group":"get_started","trigger":"finished","log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"cust_acq_web_002":{"experiment_id":"7123657947794","type":"visitor","group":"control","trigger":"hash_visitor","schedule_ts":1716500056,"log_exposures":true,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"activation_browser_deprecation_warning_may_2024":{"experiment_id":"7065959527591","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1715880596,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"activation_enterprise_signin_primer":{"experiment_id":"6443324713893","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1709667232,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"crypto_sidecar_for_comparable_keychains":{"experiment_id":"7041197731428","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1715184970,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"features_meta_rd2":{"experiment_id":"7088377708401","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1715114558,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"new_paid_lp":{"experiment_id":"6818768684695","type":"visitor","group":"treatment","trigger":"finished","log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"slack_elevate_launch":{"experiment_id":"6966627699558","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1713959798,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"cust_acq_slack_ai_apr_17":{"experiment_id":"6946248676101","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1713441387,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"marketing_recaptcha_hc":{"experiment_id":"6963734115829","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1713301140,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"marketing_hc_flow_specifier":{"experiment_id":"6989238991504","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1713296815,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"app_directory_connectors_collection":{"experiment_id":"6321714753558","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1705448247,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"anthony_test_visitor_1":{"experiment_id":"6823470010164","type":"visitor","group":"treatment","trigger":"launch_visitor","schedule_ts":1712613615,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"contact_sales_dept_removal":{"experiment_id":"6538486873169","type":"visitor","group":"treatment","trigger":"hash_visitor","schedule_ts":1711997993,"log_exposures":true,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"slack_developer_program_tos":{"experiment_id":"6725223824534","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1710376976,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"expanded_footer":{"experiment_id":"6586385434133","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1710262575,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"desktop_updater_v3_public_latest_release":{"experiment_id":"6724457596097","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1709680952,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"sk_compliance_24":{"experiment_id":"6711839338566","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1709761959,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"skip_supercookie_rotate_when_removed":{"experiment_id":"5239448035205","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1684284891,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"verify_session_from_supercookie":{"experiment_id":"5335820669477","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1686073418,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"slack_developer_program":{"experiment_id":"5782848233798","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1709242483,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"marketing_media_kit":{"experiment_id":"6696687337684","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1709232747,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"clean_up_users_and_users_team_table":{"experiment_id":"5443336696086","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1692367247,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"valid_global_session_when_fetch_cookie_sessions":{"experiment_id":"5771070997317","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1695047294,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"on24_extension":{"experiment_id":"4772019824211","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1675891595,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"slack_phone_number_exp":{"experiment_id":"5760998465057","type":"visitor","group":"treatment","trigger":"finished","log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"canvas_templates_pricing":{"experiment_id":"6538548631473","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1707410717,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"mris_usercreator_live":{"experiment_id":"4921780175444","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1706216435,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"mris_extension":{"experiment_id":"4746206947365","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1706216371,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"id_id_policy":{"experiment_id":"6488828265426","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1706045086,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"workato_form_enterprise_demo":{"experiment_id":"6457731483265","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1705509651,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"workato_form_certified_demo":{"experiment_id":"6446734501494","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1705509622,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"app_directory_connectors":{"experiment_id":"6144504493874","type":"visitor","group":"treatment","trigger":"launch_visitor","schedule_ts":1705354312,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"out_of_office_xmas_jp":{"experiment_id":"6296845198293","type":"visitor","group":"off","trigger":"finished","log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"out_of_office_xmas":{"experiment_id":"6322553087328","type":"visitor","group":"off","trigger":"finished","log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"marketing_hreflang_errors_fix":{"experiment_id":"6319747700807","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1702931766,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"eg_pricing":{"experiment_id":"6266727458225","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1702587412,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"new_gated_demo":{"experiment_id":"6171698537921","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1701209093,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"ia4_asset_refresh":{"experiment_id":"5640552621317","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1700246522,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"ia4_t3_asset_refresh":{"experiment_id":"6177775777761","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1700246503,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"marketing_cj":{"experiment_id":"5820701519667","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1699033035,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"deny_russian_ip":{"experiment_id":"3201051153989","type":"visitor","group":"on","trigger":"finished","log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"pxp862visitor":{"experiment_id":"5396563856918","type":"visitor","group":"on","trigger":"finished","log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"asean_website_personalization":{"experiment_id":"5937603225537","type":"visitor","group":"treatment","trigger":"launch_visitor","schedule_ts":1696267384,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"activation_browser_deprecation_red_banner_september_2023":{"experiment_id":"5027493502468","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1691614144,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"trials_tabs":{"experiment_id":"5755848165412","type":"visitor","group":"control","trigger":"finished","log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"expanded_nav_tab_mob":{"experiment_id":"5454526065121","type":"visitor","group":"treatment","trigger":"launch_visitor","schedule_ts":1690319627,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"expanded_nav_desktop":{"experiment_id":"5311612083732","type":"visitor","group":"treatment","trigger":"launch_visitor","schedule_ts":1684852542,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"swap_ukraine_logo_toggle":{"experiment_id":"5598910456034","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1689885040,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"contact_sales_ux":{"experiment_id":"4980218157648","type":"visitor","group":"control","trigger":"finished","log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"customer_awards_launch":{"experiment_id":"2673548411155","type":"visitor","group":"on","trigger":"finished","log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"slack_ie_address":{"experiment_id":"4857849748754","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1677793396,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"proj_brand_customer_stories_lp":{"experiment_id":"3448021380448","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1653596127,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"govslack_provision_block":{"experiment_id":"3289482788501","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1650990072,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"search_zd_vs_solr":{"experiment_id":"1355709002145","type":"visitor","group":"control","trigger":"finished","log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"digital_first_lightning_strike_custacq":{"experiment_id":"2220660679364","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1625075563,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"cust_acq_partners_template":{"experiment_id":"2232204551504","type":"visitor","group":"treatment","trigger":"launch_visitor","schedule_ts":1628191410,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"community_launch":{"experiment_id":"2652841576373","type":"visitor","group":"on","trigger":"launch_visitor","schedule_ts":1635871147,"log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"app_directory_partner_promotion_banners":{"experiment_id":"3009458145893","type":"visitor","group":"control","trigger":"finished","log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"},"app_directory_aws_promotion_banner":{"experiment_id":"3025397781073","type":"visitor","group":"control","trigger":"finished","log_exposures":false,"exposure_id":"186cc45b3cee1d025a755beeb14bca51"}},"no_login":false};</script><script type="text/javascript" crossorigin="anonymous" src="https://a.slack-edge.com/bv1-13/primer-vendor.bc1cf70043ab8c1aac02.primer.min.js" onload="window._cdn ? _cdn.ok(this, arguments) : null" onerror="window._cdn ? _cdn.failed(this, arguments) : null"></script><script type="text/javascript" crossorigin="anonymous" src="https://a.slack-edge.com/bv1-13/login-core.d683b8d740d514f731d7.primer.min.js" onload="window._cdn ? _cdn.ok(this, arguments) : null" onerror="window._cdn ? _cdn.failed(this, arguments) : null"></script><link href="https://a.slack-edge.com/bv1-13/login-core.e7910a12c731f2b97291.primer.min.css" rel="stylesheet" type="text/css" onload="window._cdn ? _cdn.ok(this, arguments) : null" onerror="window._cdn ? _cdn.failed(this, arguments) : null" crossorigin="anonymous"><link href="https://a.slack-edge.com/1b7acd1/style/rollup-slack_kit_base.css" rel="stylesheet" id="slack_kit_helpers" type="text/css" onload="window._cdn ? _cdn.ok(this, arguments) : null" onerror="window._cdn ? _cdn.failed(this, arguments) : null" crossorigin="anonymous"><link href="https://a.slack-edge.com/1789af2/style/rollup-slack_kit_helpers.css" rel="stylesheet" id="slack_kit_helpers" type="text/css" onload="window._cdn ? _cdn.ok(this, arguments) : null" onerror="window._cdn ? _cdn.failed(this, arguments) : null" crossorigin="anonymous">

<!-- slack-www-hhvm-main-iad-fgqx/ 2024-05-27 23:47:22/ vd3fa707011fbfd2dbb2c2aad4e49df01a3e4d4fd/ B:H -->

</body></html>