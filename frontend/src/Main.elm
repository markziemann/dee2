module Main exposing (..)

import Browser exposing (Document)
import Browser.Navigation as Nav
import Html exposing (..)
import Html.Attributes exposing (..)
import Info exposing (introduction)
import Maybe.Extra
import Nav exposing (navbar)
import ResultsPage.Helpers exposing (updateSearchData)
import ResultsPage.Main as RPMain
import ResultsPage.Views exposing (viewSearchResults)
import Routes
import SearchPage.Helpers exposing (defaultSearchParameters, withPagination)
import SearchPage.Main as SPMain
import SearchPage.Types exposing (SearchParameters(..))
import SearchPage.Views exposing (viewLargeSearchBar, viewSearchButton, viewSearchModeSelector)
import Types exposing (..)
import Url


initializeModelTORoute : Model -> ( Model, Cmd Msg )
initializeModelTORoute model =
    case model.route of
        Routes.SearchRoute searchUrlParameters ->
            update (RequestSearch searchUrlParameters) model

        _ ->
            ( model, Cmd.none )


init : () -> Url.Url -> Nav.Key -> ( Model, Cmd Msg )
init flags url navKey =
    let
        route =
            Routes.determinePage url

        model =
            { navKey = navKey
            , url = url
            , searchPage = SPMain.init
            , resultsPage = RPMain.init
            , route = route
            }
    in
    initializeModelTORoute model



---- UPDATE ----


unwrapMaybeMap : (a -> b -> a) -> a -> Maybe b -> a
unwrapMaybeMap function a maybeB =
    case maybeB of
        Just b ->
            function a b

        Nothing ->
            a


update : Msg -> Model -> ( Model, Cmd Msg )
update msg model =
    let
        fromSearchPage =
            \( mdl, cmd, maybeOutMsg ) ->
                ( { model
                    | searchPage = mdl
                    , resultsPage = unwrapMaybeMap updateSearchData model.resultsPage maybeOutMsg
                  }
                , Maybe.Extra.values
                    [ Just <| Cmd.map GotSearchPageMsg cmd
                    , Maybe.map
                        (.searchParameters
                            >> Routes.searchResultsRoute
                            >> Nav.pushUrl model.navKey
                        )
                        maybeOutMsg
                    ]
                    |> Cmd.batch
                )

        fromResultsPage =
            \( mdl, cmd, maybePaginationOffset ) ->
                ( { model
                    | resultsPage = mdl
                  }
                , Maybe.Extra.values
                    [ Just <| Cmd.map GotResultsPageMsg cmd
                    , Maybe.map
                        ((withPagination <| defaultSearchParameters model.searchPage)
                            >> Routes.searchResultsRoute
                            >> Nav.pushUrl model.navKey
                        )
                        maybePaginationOffset
                    ]
                    |> Cmd.batch
                )
    in
    case msg of
        GotSearchPageMsg message ->
            SPMain.update message model.searchPage
                |> fromSearchPage

        GotResultsPageMsg message ->
            RPMain.update message model.resultsPage
                |> fromResultsPage

        RequestSearch searchParameters ->
            SPMain.update
                (SPMain.search searchParameters)
                model.searchPage
                |> fromSearchPage

        LinkClicked urlRequest ->
            case urlRequest of
                Browser.Internal url ->
                    ( model, Nav.pushUrl model.navKey (Url.toString url) )

                Browser.External href ->
                    ( model, Nav.load href )

        HomeReset ->
            ( { model | searchPage = SPMain.init }, Nav.pushUrl model.navKey "/" )

        UrlChanged url ->
            ( { model | url = url, route = Routes.determinePage url }, Cmd.none )



---- VIEW ----


pageLayout : List (Html Msg) -> List (Html Msg)
pageLayout content =
    [ navbar
    , div [ class "container my-5 mx-auto" ] content
    , introduction
    ]


pageView : Model -> List (Html Msg)
pageView model =
    let
        fromSearchPage =
            List.map (Html.map GotSearchPageMsg)

        fromResultsPage =
            List.map (Html.map GotResultsPageMsg)
    in
    pageLayout <|
        case model.route of
            Routes.HomeRoute ->
                fromSearchPage
                    [ viewLargeSearchBar model.searchPage
                    , viewSearchModeSelector model.searchPage.searchMode
                    , viewSearchButton model.searchPage
                    ]

            Routes.SearchRoute (SearchParameters _ _ paginationOffset) ->
                fromResultsPage
                    (viewSearchResults model.resultsPage paginationOffset)

            Routes.Unknown ->
                [ text "Hmm... I don't recognise that url." ]


view : Model -> Document Msg
view model =
    { title = "Digital Expression Explorer 2"
    , body = pageView model
    }


subscriptions : Model -> Sub Msg
subscriptions model =
    case model.route of
        Routes.HomeRoute ->
            Sub.batch
                [ Sub.map GotSearchPageMsg <| SPMain.subscriptions model.searchPage
                ]

        Routes.SearchRoute _ ->
            Sub.none

        Routes.Unknown ->
            subscriptions { model | route = Routes.HomeRoute }



---- PROGRAM ----


main : Program () Model Msg
main =
    Browser.application
        { view = view
        , init = init
        , update = update
        , subscriptions = subscriptions
        , onUrlRequest = LinkClicked
        , onUrlChange = UrlChanged
        }
