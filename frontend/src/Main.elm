module Main exposing (..)

import Browser exposing (Document)
import Browser.Navigation as Nav
import Html exposing (..)
import Html.Attributes exposing (..)
import Info exposing (introduction)
import Maybe.Extra
import Nav exposing (navbar)
import ResultsPage.Main as RPMain
import ResultsPage.Types
import ResultsPage.Views exposing (viewSearchResults)
import Routes
import SearchPage.Helpers exposing (getSearchMode, sameModeAndString, withPagination)
import SearchPage.Main as SPMain
import SearchPage.Types exposing (SearchParameters(..))
import SearchPage.Views exposing (viewLargeSearchBar, viewSearchButton, viewSearchModeSelector)
import SharedTypes exposing (PaginationOffset, RemoteData(..))
import Types exposing (..)
import Url


initializeModelTORoute : Model -> ( Model, Cmd Msg )
initializeModelTORoute model =
    case model.route of
        Routes.ResultsRoute searchUrlParameters ->
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
            , resultsPage = RPMain.init NotAsked Nothing
            , route = route
            }
    in
    initializeModelTORoute model



---- UPDATE ----


maybeUpdateResults : ResultsPage.Types.Model -> Maybe SearchPage.Types.OutMsg -> ResultsPage.Types.Model
maybeUpdateResults model maybeOutMsg =
    case maybeOutMsg of
        Just { searchResults, searchParameters } ->
            case model.searchParameters of
                Just currentSearchParams ->
                    if sameModeAndString searchParameters currentSearchParams then
                        { model | searchResults = searchResults, searchParameters = Just searchParameters }

                    else
                        RPMain.init searchResults (Just searchParameters)

                _ ->
                    RPMain.init searchResults (Just searchParameters)

        _ ->
            model


maybeUpdateUrl cmd model =
    .searchParameters
        >> Routes.searchResultsRoute
        >> Nav.pushUrl model.navKey


update : Msg -> Model -> ( Model, Cmd Msg )
update msg model =
    let
        fromSearchPage =
            \( mdl, cmd, maybeOutMsg ) ->
                ( { model
                    | searchPage = mdl
                    , resultsPage = maybeUpdateResults model.resultsPage maybeOutMsg
                  }
                , case maybeOutMsg of
                    Just { searchResults, searchParameters } ->
                        if searchResults == Loading then
                            Debug.log "PUSH" Cmd.batch
                                [ Cmd.map GotSearchPageMsg cmd
                                , Nav.pushUrl model.navKey <| Routes.searchResultsRoute searchParameters
                                ]

                        else
                            Cmd.map GotSearchPageMsg cmd

                    Nothing ->
                        Cmd.map GotSearchPageMsg cmd
                )

        fromResultsPage =
            \( mdl, cmd, maybePaginationOffset ) ->
                let
                    newModel =
                        { model | resultsPage = mdl }

                    command =
                        Cmd.map GotResultsPageMsg cmd
                in
                case ( maybePaginationOffset, mdl.searchParameters ) of
                    ( Just paginationOffset, Just searchParameters ) ->
                        update (RequestSearch <| withPagination paginationOffset searchParameters) newModel
                            |> (\( m, c ) -> ( m, Cmd.batch [ c, command ] ))

                    ( _, _ ) ->
                        ( newModel, command )
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
                    Debug.log "Second" ( model, Nav.pushUrl model.navKey (Url.toString url) )

                Browser.External href ->
                    ( model, Nav.load href )

        HomeReset ->
            ( { model | searchPage = SPMain.init }, Nav.pushUrl model.navKey "/" )

        UrlChanged url ->
            let
                route =
                    Routes.determinePage url

                default =
                    ( { model | url = url, route = route }, Cmd.none )
            in
            case ( route, model.resultsPage.searchParameters ) of
                ( Routes.ResultsRoute searchParameters, Nothing ) ->
                    -- This can occur when the browser navigates forwards after
                    -- Leaving the web site
                    update (RequestSearch searchParameters) model

                ( Routes.ResultsRoute searchParameters, Just searchParams ) ->
                    if searchParameters /= searchParams then
                        update (RequestSearch searchParameters) model

                    else
                        default

                ( _, _ ) ->
                    default



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
                    , viewSearchModeSelector <| getSearchMode model.searchPage.searchParameters
                    , viewSearchButton model.searchPage
                    ]

            Routes.ResultsRoute (SearchParameters _ _ paginationOffset) ->
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

        Routes.ResultsRoute _ ->
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
